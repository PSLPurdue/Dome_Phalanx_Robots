function [UTOT, grad_E] = sym_elasticenergy(dof, Ks, Bis_S, coords0, connecRigid, connecBistable, connecTorsion, BC_nodes, Nnodes)
% SYM_ELASTICENERGY   Compute total elastic energy and its gradient
%   Inputs:
%     dof            – vector of free-node displacements
%     Ks             – [M×3] array: [k_rigid, k_lim, ktorsion]
%     Bis_S          – [M×3] array: [kBS, al, dd]
%     coords0        – [2×Nnodes] original node coordinates
%     connecRigid    – [2P×2] rigid & limit-layer spring pairs
%     connecBistable – [Q×2] bistable spring pairs
%     connecTorsion  – [R×3] torsional spring triplets
%     BC_nodes       – indices of fixed nodes
%     Nnodes         – total number of nodes
%   Outputs:
%     UTOT           – scalar total elastic energy
%     grad_E         – (optional) row vector gradient of energy

%% 1. Define energy functions for each spring type
% 1.1 Bistable spring energy (double-well form)
Ubs = @(xij, x0ij, kBS, al, dd) 0.5 .* kBS .* (xij - x0ij).^2 .* ...
    (1 + (1 - al) .* ( ((xij - x0ij)./dd).^2 - 2 .* ((xij - x0ij)./dd) ));
% 1.2 Linear spring energy (rigid & limit-layer)
UL = @(xij, x0ij, k) 0.5 * k .* (xij - x0ij).^2;
% 1.3 Torsional spring energy
UTorsion = @(theta, theta0, kt) 0.5 * kt .* (theta - theta0).^2;

%% 2. Map DOF vector back to full displacement matrix
% 2.1 Identify free nodes
NoBcNodes = 1:Nnodes;
NoBcNodes(BC_nodes) = [];
% 2.2 Reshape displacement vector to [2×#freeNodes]
uvw = zeros(2, Nnodes);
uvw(:, NoBcNodes) = reshape(dof, [2, length(NoBcNodes)]);

%% 3. Unpack spring constants and bistable parameters
k_rigid   = Ks(:,1);          % rigid spring stiffness
k_lim     = Ks(:,2);          % limit-layer spring stiffness
ktorsion  = Ks(:,3);          % torsional spring stiffness
kBS       = Bis_S(:,1);       % bistable spring stiffness
al        = Bis_S(:,2);       % bistable shape parameter
dd        = Bis_S(:,3);       % displacement scale for bistability
% 3.1 Build array of torsion constants for each hinge edge
temp = reshape([ktorsion; ktorsion], 1, []);
ktorsion_array = temp(1:end-1);
% 3.2 Build array for alternating rigid springs
k_rigid_array = [k_rigid; k_rigid(1)];

%% 4. Compute deformed coordinates
coords = coords0 + uvw;  % add displacements to original positions

%% 5. Compute distances for rigid springs
% 5.1 Extract indices for rigid springs (odd rows)
connecRigidi = connecRigid(1:2:end, 1);
connecRigidj = connecRigid(1:2:end, 2);
% 5.2 Rest lengths
XD0 = coords0(:, connecRigidi) - coords0(:, connecRigidj);
lijR0 = sqrt(sum(XD0.^2))';
% 5.3 Current lengths
XD = coords(:, connecRigidi) - coords(:, connecRigidj);
lijR = sqrt(sum(XD.^2))';

%% 6. Compute distances for limit-layer springs
% 6.1 Extract indices (even rows)
connecLimi = connecRigid(2:2:end, 1);
connecLimj = connecRigid(2:2:end, 2);
% 6.2 Rest lengths
XD0 = coords0(:, connecLimi) - coords0(:, connecLimj);
lijlim0 = sqrt(sum(XD0.^2))';
% 6.3 Current lengths
XD = coords(:, connecLimi) - coords(:, connecLimj);
lijlim = sqrt(sum(XD.^2))';

%% 7. Compute distances for bistable springs
connecBistablei = connecBistable(:,1);
connecBistablej = connecBistable(:,2);
% Rest lengths
XD0 = coords0(:, connecBistablei) - coords0(:, connecBistablej);
lijB0 = sqrt(sum(XD0.^2))';
% Current lengths
XD = coords(:, connecBistablei) - coords(:, connecBistablej);
lijB = sqrt(sum(XD.^2))';

%% 8. Compute angles for torsional springs
nTsprings = size(connecTorsion, 1);
theta0 = zeros(nTsprings, 1);
theta  = zeros(nTsprings, 1);
for n = 1:nTsprings
    % 8.1 Reference positions
    x_i_0 = coords0(:, connecTorsion(n,1))';
    x_j_0 = coords0(:, connecTorsion(n,2))';
    x_k_0 = coords0(:, connecTorsion(n,3))';
    % Edge vectors in rest config
    e0_0 = x_j_0 - x_i_0;
    e1_0 = x_k_0 - x_i_0;
    % Normalize
    e0_0_h = e0_0 / norm(e0_0);
    e1_0_h = e1_0 / norm(e1_0);
    % Rest angle
    theta0(n) = acos(dot(e0_0_h, e1_0_h));
    
    % 8.2 Current positions
    x_i = coords(:, connecTorsion(n,1))';
    x_j = coords(:, connecTorsion(n,2))';
    x_k = coords(:, connecTorsion(n,3))';
    % Edge vectors
    e0 = x_j - x_i;
    e1 = x_k - x_i;
    % Normalize
    e0_h = e0 / norm(e0);
    e1_h = e1 / norm(e1);
    % Current angle
    theta(n) = acos(dot(e0_h, e1_h));
end

%% 9. Sum energy contributions
% 9.1 Rigid springs energy
UR   = sum(UL(lijR,   lijR0,   k_rigid_array));
% 9.2 Limit-layer energy
Ulim = sum(UL(lijlim, lijlim0, k_lim));
% 9.3 Bistable springs energy
UBS  = sum(Ubs(lijB,   lijB0,   kBS, al, dd));
% 9.4 Torsional springs energy
UT   = sum(UTorsion(theta', theta0', ktorsion_array));
% 9.5 Total energy
UTOT = UR + Ulim + UBS + UT;

%% 10. Compute gradient if requested
if nargout > 1
    grad_E = Grad_sym_elasticenergy(dof, Ks, Bis_S, coords0, ...
        connecRigid, connecBistable, connecTorsion, BC_nodes, Nnodes);
end
end

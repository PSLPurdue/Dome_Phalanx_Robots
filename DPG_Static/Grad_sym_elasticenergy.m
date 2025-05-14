function grad_E = Grad_sym_elasticenergy(dof, Ks, Bis_S, coords0, connecRigid, connecBistable, connecTorsion, BC_nodes, Nnodes)
% GRAD_SYM_ELASTICENERGY  Compute gradient of total elastic energy
%   Given degrees of freedom and spring definitions, assemble internal forces
%   from rigid, limit-layer, bistable, and torsional springs, then return
%   the force vector for free nodes (negative energy gradient).
%
%   Inputs:
%     dof            – [2*(Nnodes-#BC)×1] vector of free-node displacements
%     Ks             – [M×3] array of linear spring constants [k_rigid, k_lim, ktorsion]
%     Bis_S          – [M×3] array of bistable params [kBS, alpha, dd]
%     coords0        – [2×Nnodes] matrix of undeformed node coordinates
%     connecRigid    – [2P×2] rigid & limit-layer spring node pairs
%     connecBistable – [Q×2] bistable spring node pairs
%     connecTorsion  – [R×3] triples defining torsional springs (vertex, neighbor1, neighbor2)
%     BC_nodes       – indices of nodes with fixed boundary conditions
%     Nnodes         – total number of nodes in the system
%
%   Output:
%     grad_E         – [1×2*(Nnodes-#BC)] row vector of internal forces for free DOFs

%% 1. Map DOFs back to full coordinate array
% Determine free-node indices
allNodes    = 1:Nnodes;
freeNodes   = setdiff(allNodes, BC_nodes);
% Reshape dof vector into [2×#freeNodes] displacement matrix
uvw = zeros(2, Nnodes);
uvw(:, freeNodes) = reshape(dof, [2, numel(freeNodes)]);

%% 2. Unpack spring constants and bistable params
k_rigid    = Ks(:,1);           % stiff rigid springs
k_lim      = Ks(:,2);           % limit-layer springs
ktorsion   = Ks(:,3);           % torsional spring constants
kBS        = Bis_S(:,1);        % bistable spring stiffness
alpha      = Bis_S(:,2);        % bistable shape parameter
dd         = Bis_S(:,3);        % displacement scale for bistability
% Build arrays for torsion along each edge
kt_arr = repelem(ktorsion, 1, 2);
% Wrap rigid-array for alternating connections
k_rig_arr = [k_rigid; k_rigid(1:end)];

%% 3. Compute deformed node positions
coords = coords0 + uvw;  % add displacement to original coords

%% 4. Initialize internal force accumulator
F_in = zeros(size(coords));  % 2×Nnodes

%% 5. Rigid-body springs (even rows in connecRigid)
% Extract node indices for rigid springs
iR = connecRigid(1:2:end,1);
jR = connecRigid(1:2:end,2);
% Rest lengths
XR0 = coords0(:, iR) - coords0(:, jR);
s0R = sqrt(sum(XR0.^2,1));
for idx = 1:numel(iR)
    % current separation
    dX = coords(:, iR(idx)) - coords(:, jR(idx));
    curL = norm(dX);
    % linear spring force = k*(1 - rest/cur)*dX
    F = k_rig_arr(idx) * (1 - s0R(idx)/curL) * dX;
    % distribute to nodes
    F_in(:, iR(idx)) = F_in(:, iR(idx)) + F;
    F_in(:, jR(idx)) = F_in(:, jR(idx)) - F;
end

%% 6. Limit-layer springs (odd rows in connecRigid)
iL = connecRigid(2:2:end,1);
jL = connecRigid(2:2:end,2);
XL0 = coords0(:, iL) - coords0(:, jL);
s0L = sqrt(sum(XL0.^2,1));
for idx = 1:numel(iL)
    dX = coords(:, iL(idx)) - coords(:, jL(idx));
    curL = norm(dX);
    F = k_lim(idx) * (1 - s0L(idx)/curL) * dX;
    F_in(:, iL(idx)) = F_in(:, iL(idx)) + F;
    F_in(:, jL(idx)) = F_in(:, jL(idx)) - F;
end

%% 7. Bistable springs
iB = connecBistable(:,1);
jB = connecBistable(:,2);
XB0 = coords0(:, iB) - coords0(:, jB);
s0B = sqrt(sum(XB0.^2,1));
for idx = 1:numel(iB)
    dX = coords(:, iB(idx)) - coords(:, jB(idx));
    ext = norm(dX) - s0B(idx);  % extension beyond rest
    % non-linear bistable force
    Fmag = kBS(idx) * ext * (1 + (1-alpha(idx))*(2*(ext/dd(idx))^2 - 3*(ext/dd(idx))));
    dir = dX / norm(dX);
    % distribute
    F = Fmag * dir;
    F_in(:, iB(idx)) = F_in(:, iB(idx)) + F;
    F_in(:, jB(idx)) = F_in(:, jB(idx)) - F;
end

%% 8. Torsional springs (angle at node i between j and k)
nT = size(connecTorsion,1);
for idx = 1:nT
    % reference triangle
    tri = connecTorsion(idx,:);
    x0 = coords0(:,tri(1))'; x1 = coords0(:,tri(2))'; x2 = coords0(:,tri(3))';
    % edges in reference
    e0r = x1 - x0; e1r = x2 - x0;
    % rest-angle
    th0 = acos(dot(e0r/norm(e0r), e1r/norm(e1r)));
    % current triangle
    y0 = coords(:,tri(1))'; y1 = coords(:,tri(2))'; y2 = coords(:,tri(3))';
    e0 = y1 - y0; e1 = y2 - y0;
    a = norm(e0); b = norm(e1);
    % current angle
    th = acos(dot(e0/norm(e0), e1/norm(e1)));
    % dot product
    de = dot(e0, e1);
    % gradients per formula
    grad0 = -(e0+e1)/(a*b) + de*(e0/(a^3*b) + e1/(a*b^3));
    grad1 =  e1/(a*b) - (de/(a^3*b))*e0;
    grad2 =  e0/(a*b) - (de/(a*b^3))*e1;
    % torsional moment → force
    factor = -kt_arr(idx) * ((th-th0)/sin(th));
    F_in(:,tri(1)) = F_in(:,tri(1)) + factor*grad0';
    F_in(:,tri(2)) = F_in(:,tri(2)) + factor*grad1';
    F_in(:,tri(3)) = F_in(:,tri(3)) + factor*grad2';
end

%% 9. Extract forces for free nodes and reshape
F_free = F_in(:, freeNodes);
grad_E = reshape(F_free, 1, []);  % row vector
end
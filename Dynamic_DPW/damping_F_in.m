function TOT = damping_F_in(dof, dof_p, eta, coords0, connecBistable, BC_nodes, Nnodes)
% DAMPING_F_IN  Compute damping forces for bistable springs
%   Constructs viscous damping force vector based on relative node velocities
%
%   Inputs:
%     dof            – displacement vector for free DOFs
%     dof_p          – velocity vector for free DOFs
%     eta            – damping coefficient
%     coords0        – [2×Nnodes] original node coordinates
%     connecBistable – [Q×2] pairs of nodes connected by bistable springs
%     BC_nodes       – indices of fixed (boundary) nodes
%     Nnodes         – total number of nodes
%
%   Output:
%     TOT            – [2*(Nnodes-#BC)×1] damping force vector for free DOFs

%% 1. Initialize arrays for displacement, velocity, and damping force
uvw        = zeros(size(coords0));   % displacement field (2×Nnodes)
vel        = zeros(size(coords0));   % velocity field (2×Nnodes)
damp_force = zeros(size(coords0));   % accumulated damping forces

%% 2. Map DOF vectors to full nodal fields
% 2.1 Identify free (non-BC) node indices
NoBcNodes = 1:Nnodes;
NoBcNodes(BC_nodes) = [];
% 2.2 Reshape displacement and velocity into [2×#freeNodes]
uvw(:, NoBcNodes) = reshape(dof,   [2, length(NoBcNodes)]);
vel(:, NoBcNodes) = reshape(dof_p, [2, length(NoBcNodes)]);

%% 3. Compute current deformed coordinates
coords = coords0 + uvw;  % deformed positions (2×Nnodes)

%% 4. Precompute rest and current distances for bistable springs
% 4.1 Node pairs for bistable springs
iB = connecBistable(:,1);
jB = connecBistable(:,2);
% 4.2 Rest distances (sij)
X0 = coords0(:, iB) - coords0(:, jB);
sij = sqrt(sum(X0.^2))';
% 4.3 Current distances
x_ij = coords(:, iB) - coords(:, jB);
norm_x_ij = sqrt(sum(x_ij.^2))';

%% 5. Gather nodal velocities for each spring endpoint
Vi = vel(:, iB);  % velocities at node i of each spring
Vj = vel(:, jB);  % velocities at node j of each spring

%% 6. Compute damping force contribution for each spring
% Formula: f_d = eta*[ (1 - sij/norm_x)*relative_vel + (sij/norm_x^3)*(x·relative_vel)*x ]
for k = 1:length(iB)
    rel_vel = Vi(:,k) - Vj(:,k);         % relative velocity vector
    x_vec   = x_ij(:,k);                % current spring vector
    x_norm  = norm_x_ij(k);
    s0      = sij(k);
    % viscous damping force magnitude & direction
    f_d_ij = eta * ( ...
        (1 - s0/x_norm) * rel_vel + ...                     % tangential part
        (s0/(x_norm^3)) * (x_vec' * rel_vel) * x_vec ...   % radial part
    );
    % distribute half to each node
    damp_force(:, iB(k)) = damp_force(:, iB(k)) + f_d_ij;
    damp_force(:, jB(k)) = damp_force(:, jB(k)) - f_d_ij;
end

%% 7. Extract damping forces for free DOFs and reshape
F_free = damp_force(:, NoBcNodes);  % forces only at free nodes
numVar = length(NoBcNodes) * 2;      % total free DOFs
TOT    = reshape(F_free, [numVar, 1]);  % column vector for solver

end
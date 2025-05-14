function F_BC = F_ext(t, x, P0, tau_1, tau_2, tau_3, tau_4, k_dome, connecRigid, connecBistable, BC_nodes, coords0, Nnodes)
% F_EXT  Compute external force vector for time-varying pressure loads
%   Applies normal forces on limiting, bistable, and rigid layers
%   based on time-dependent pressure P_t, then returns forces on free DOFs.
%
%   Inputs:
%     t               – current time
%     x               – solution vector (dofs) for all free nodes
%     P0, tau_1-4     – parameters for time-dependent pressure P_t
%     k_dome          – dome stiffness multiplier for rigid springs
%     connecRigid     – [2P×2] rigid & limiting layer connectivity
%     connecBistable  – [Q×2] bistable spring connectivity
%     BC_nodes        – indices of fixed (boundary) nodes
%     coords0         – [2×Nnodes] original undeformed coordinates
%     Nnodes          – total number of nodes
%
%   Output:
%     F_BC            – [2*(Nnodes-#BC)×1] external force vector for free DOFs

%% 1. Initialize displacement & force arrays
uvw = zeros(size(coords0));    % full displacement field (2×Nnodes)
F   = zeros(size(coords0));    % accumulated force field (2×Nnodes)

% 1.1 Compute current pressure from prescribed profile
P_time_i = P_t(P0, t, tau_1, tau_2, tau_3, tau_4);

% 1.2 Map solution vector x to full node displacements
NoBcNodes = 1:Nnodes;
NoBcNodes(BC_nodes) = [];               % indices of free nodes
uvw(:, NoBcNodes) = reshape(x, [2, length(NoBcNodes)]);
coords_deformed = coords0 + uvw;        % deformed coordinates (2×Nnodes)

%% 2. Identify element groups (layers)
lim_layer    = connecRigid(2:2:end, :); % limiting layer springs (even rows)
dome_layer   = connecBistable;         % bistable dome springs
rigid_connect= connecRigid(1:2:end, :); % rigid springs (odd rows)

%% 3. Limiting layer: compute normal forces
normal_vec_lim = zeros(size(lim_layer,1), 2);
for i = 1:size(lim_layer,1)
    % 3.1 Compute tangent vector and face area A_layer
    tangent = coords_deformed(:, lim_layer(i,2)) - coords_deformed(:, lim_layer(i,1));
    A_layer = norm(tangent);
    % 3.2 Compute outward normal (rotate tangent by +90°)
    n = [-tangent(2), tangent(1)];
    n_h = -n / norm(n);  % unit normal (negated for consistent orientation)
    % 3.3 Scale by pressure and face area
    normal_vec_lim(i, :) = P_time_i * A_layer * n_h;
end
% 3.4 Distribute half force to each node of the face
F(:, lim_layer(:,1)) = F(:, lim_layer(:,1)) + normal_vec_lim'/2;
F(:, lim_layer(:,2)) = F(:, lim_layer(:,2)) + normal_vec_lim'/2;

%% 4. Bistable (dome) springs: compute normal forces
normal_vec_dome = zeros(size(dome_layer,1), 2);
for i = 1:size(dome_layer,1)
    tangent = coords_deformed(:, dome_layer(i,2)) - coords_deformed(:, dome_layer(i,1));
    A_dome  = norm(tangent);
    n       = [-tangent(2), tangent(1)];
    n_h     = n / norm(n);  % unit normal (no negation)
    normal_vec_dome(i, :) = P_time_i * A_dome * n_h;
end
F(:, dome_layer(:,1)) = F(:, dome_layer(:,1)) + normal_vec_dome'/2;
F(:, dome_layer(:,2)) = F(:, dome_layer(:,2)) + normal_vec_dome'/2;

%% 5. Rigid springs: compute normal forces (skip first spring)
normal_vec_rigid = zeros(size(rigid_connect,1), 2);
for i = 2:size(rigid_connect,1)
    tangent = coords_deformed(:, rigid_connect(i,2)) - coords_deformed(:, rigid_connect(i,1));
    A_rigid = norm(tangent);
    n       = [-tangent(2), tangent(1)];
    n_h     = -n / norm(n);  % unit normal
    normal_vec_rigid(i, :) = k_dome * P_time_i * A_rigid * n_h;
end
F(:, rigid_connect(:,1)) = F(:, rigid_connect(:,1)) + normal_vec_rigid'/2;
F(:, rigid_connect(:,2)) = F(:, rigid_connect(:,2)) + normal_vec_rigid'/2;

%% 6. Extract forces for free nodes and reshape for solver
F_BC = F(:, NoBcNodes);  % keep only columns for free DOFs
numVar = length(NoBcNodes)*2;
F_BC = reshape(F_BC, [numVar, 1]);  % column vector for solver

end

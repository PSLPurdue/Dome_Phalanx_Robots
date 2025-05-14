function [time, x, connecNodes, Nnodes, BC_nodes, coords0] = section_i_robot(Hv, Ch_L, Unit_sep, t, t_lim, rb, which_state, P0, tau_1, tau_2, tau_3, tau_4, tspan)
% SECTION_I_ROBOT  Setup and run dynamic simulation of DPG robot
%   This function constructs geometry, material, and dynamic parameters,
%   computes the initial inverted state via optimization, and then
%   integrates the equations of motion over tspan.
%
%   Inputs:
%     Hv            – [1×N] dome heights
%     Ch_L, Unit_sep, t, t_lim, rb – geometric parameters
%     which_state   – index for initial inverted configuration
%     P0, tau_1-4   – pressure profile parameters
%     tspan         – [t0 tf] time interval for simulation
%
%   Outputs:
%     time          – time vector from ODE solver
%     x             – state trajectories [disp; vel] over time
%     connecNodes   – connectivity for damping springs
%     Nnodes        – total number of nodes
%     BC_nodes      – indices of fixed nodes
%     coords0       – undeformed node coordinates

%% 1. Problem setup: geometric & material properties
% 1.1 Segment count based on dome vector length
number_of_segments = length(Hv);
% 1.2 Elastic properties
E  = 26.0;            % Young's modulus
nu = 0.35;            % Poisson ratio
% 1.3 Unit-cell dimensions
t_ch = 1.0;            % transition thickness
base = 2 * rb;
UC   = 1.5 * base;
% 1.4 Radius for bistability
R = (Hv.^2 + rb^2) ./ (2 * Hv);
% 1.5 Boundary conditions: first three nodes fixed
BC_nodes = [1 2 3];
% 1.6 Segment geometry
theta = 0;            % start angle
segment_length = Ch_L + t_ch + Unit_sep;
segment_height = UC/2 + t_lim + t;

%% 2. Build lattice geometry
geo_matrix = DPG_Geometry(segment_height, segment_length, number_of_segments);
coords0        = geo_matrix.coords0;         % [2×Nnodes]
connecRigid    = geo_matrix.connecRigid;     % rigid & limit-layer
connecBistable = geo_matrix.connecBistable;  % bistable springs
connecTorsion  = geo_matrix.connecTorsion;   % torsional springs
Nnodes         = geo_matrix.Nnodes;          % total nodes
connecNodes    = geo_matrix.connecNodes;     % connectivity for damping

%% 3. Compute spring constants
% 3.1 Linear springs (limiting layer)
k_lim = E * (t_lim * UC) ./ (Ch_L + t + t_ch);
k_s   = 5 * k_lim;                             % scaled stiff springs
% 3.2 Torsional springs (plate bending)
D  = (E * t_lim^3) / (12 * (1-nu^2));
kt = D * (UC ./ (Unit_sep + t + t_ch));
Ks = [k_s', k_lim', kt'];                     % [M×3]
% 3.3 Bistable spring constants
kBS = (E ./ R.^2) .* ( ...                   % nonlinear bistable stiffness
    - (1.97 .* R .* t.^5) ./ Hv.^3 + ...
    (7.4  .* Hv.^2 .* t.^3) ./ R.^2 + ...
    (3.5  .* R .* t.^4) ./ Hv.^2 + ...
    (0.37 .* Hv.^2 .* t.^2) ./ R   + ...
    (42.2 .* t.^5) ./ Hv.^2            - ...
    (35.8 .* Hv .* t.^4) ./ R.^2       + ...
    (71.8 .* t.^5) ./ (Hv .* R)        - ...
    (3.4  .* R .* t.^3) ./ Hv           - ...
    (67.7 .* t.^4) ./ Hv                + ...
    4.2  .* R .* t.^2                   + ...
    11.1 .* t.^3);

% 3.4 Shape parameter alpha (scaled by 1.3)
al = (0.3823 .* Hv.^3) ./ R.^3 - ...
     (0.5736 .* t.^3) ./ Hv.^3 - ...
     (2.8677 .* Hv.^2 .* t) ./ R.^3 - ...
     (0.9559 .* Hv.^2) ./ R.^2 - ...
     (8.9853 .* t.^3) ./ (Hv.^2 .* R) + ...
     (1.5295 .* t.^2) ./ Hv.^2 - ...
     (4.0147 .* t.^3) ./ (Hv .* R.^2) + ...
     (9.1764 .* Hv .* t) ./ R.^2 + ...
    (17.3971 .* t.^2) ./ (Hv .* R) + ...
     (0.5736 .* Hv) ./ R   + ...
    (17.7795 .* t.^3) ./ R.^3 - ...
    (19.5   .* t.^2) ./ R.^2 - ...
     (6.1177 .* t)   ./ R;

% 3.5 Tip extension offset
tip_dis = 0.2565 * t + 2.1425 * Hv;
% 3.6 Ensure non-negative dd offsets
dd = tip_dis - Hv - Unit_sep - t_ch/2 - t/2;
dd(dd < 0) = 0.01;
% 3.7 Combine bistable params
Bis_S = [kBS', al', dd'];

%% 4. Dynamic parameters
eta_iso = 0.01;         % isotropic damping
eta     = 0.05;         % viscous damping coefficient

%% 5. Mass and DOF setup
% 5.1 Degrees of freedom count
dof_nodes = Nnodes - length(BC_nodes);
numVar = 2 * dof_nodes;
% 5.2 Lumped mass per node
DPG_mass = (7e-6) * (number_of_segments/5);   % [Ton]
mass_V = (DPG_mass / Nnodes) * ones(numVar,1);

%% 6. Initial inverted state via fmincon
% 6.1 Compute approximate inverted node positions
node_locations_deformed = approx_full_inverted_state(which_state, number_of_segments, segment_length, segment_height, theta, 1.2 * dd');
node_deformations = node_locations_deformed - coords0;
deformed_coord_deltas = reshape(node_deformations(:, setdiff(1:Nnodes, BC_nodes)), [numVar,1]);
% 6.2 Optimize total elastic energy
totEn = @(d) sym_elasticenergy(d, Ks, Bis_S, coords0, connecRigid, connecBistable, connecTorsion, BC_nodes, Nnodes);
opts = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',1e6,'OptimalityTolerance',1e-6,'MaxIterations',1e5,'StepTolerance',1e-16,'SpecifyObjectiveGradient',true);
[x_sol, ~] = fmincon(totEn, deformed_coord_deltas, [],[],[],[],[],[],[], opts);

%% 7. Run dynamic simulation
% 7.1 Initial state for ODE: [optimized disp; zero vel]
x0 = [x_sol; zeros(numVar,1)];
% 7.2 Integrate using solve_system
[time, x] = solve_system(x0, tspan, numVar, mass_V, eta, eta_iso, 1.0, P0, tau_1, tau_2, tau_3, tau_4, Ks, Bis_S, connecNodes, connecRigid, connecBistable, connecTorsion, BC_nodes, coords0, Nnodes);
end

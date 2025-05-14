% Script to assemble and solve the static equilibrium of the DPG spring lattice
% Author: Juan Osorio
% Date: April 2025
%
% This script builds the geometry of a bistable spring lattice (DPG metasheet),
% defines material and geometric parameters, computes spring constants,
% then runs an optimization (fmincon) to find the equilibrium deformed state
% for objects of varying sizes. Finally, it plots the equilibrium configurations.

%% 1. Initialize Environment
clc;                % Clear command window
close all;          % Close all figure windows
clear;              % Remove all variables from workspace

%% 2. Figure Default Settings
% Set up default figure properties for publication-quality plots
set(0, ...
    'DefaultFigureColor', 'white', ...       % White background
    'DefaultAxesFontSize', 22, ...           % Larger font
    'DefaultAxesFontName', 'CMU Serif', ...  % Serif font
    'DefaultLineLineWidth', 2, ...           % Thicker lines
    'DefaultAxesBox', 'on');                 % Box around axes
% Define a color palette for plotting
colors = {'#440154', '#453781', '#2d708e', '#20a387', '#55c667', '#b8de29'};

%% 3. Object & State Definitions
% List of object sizes (e.g., heights) and corresponding precomputed states
Obj_sizeV    = [70 60 45 35];   % Four object size parameters
which_stateV = [17 25 30 31];   % Indices into a library of inverted states

% Global base orientation parameters
theta   = -27.065;              % Rotation angle [deg]
base_L  = 20.081;               % Base length

%% 4. Geometric Data (Finger / Unit Parameters)
rb     = 5;                                 % Base radius of dome
Hv     = [3.8187 4.1549 4.5246 4.8239 2.5];  % Dome heights for each unit
Ch_L   = 8.5773 * ones(size(Hv));            % Channel lengths (constant)
Unit_sep = 1.8232 * ones(size(Hv));          % Separation between units

% Thickness parameters
t     = 0.8326;                 % Dome thickness
t_lim = 1.056;                  % Limiting layer thickness

% Additional problem parameters
t_ch              = 1.0;                % Chamber thickness
number_of_segments = length(Hv);         % Number of segments in the DPF
E   = 26.0;                             % Young’s modulus [units]
nu  = 0.35;                             % Poisson’s ratio

% Derived unit-cell dimensions
base = 2 * rb;              % Base diameter
UC   = 1.5 * base;          % Unit-cell size
R    = (Hv.^2 + rb^2) ./ (2 * Hv);  % Radius of curvature of the dome

start_angle  = 0;            % Initial winding angle
BC_nodes     = [1 2 3];      % Node indices with BCs (Pin)

% Segment geometry for DPG_Geometry
segment_length = Ch_L + t_ch + Unit_sep;
segment_height = UC/2 + t_lim + t; 

%% 5. Build Lattice Geometry
% Calls external function to generate nodes & connectivity
geo_matrix      = DPG_Geometry(segment_height, segment_length, number_of_segments);
coords0         = geo_matrix.coords0;       % Original (undeformed) node coords
connecRigid     = geo_matrix.connecRigid;   % Rigid-body elements
connecBistable  = geo_matrix.connecBistable;% Bistable spring connectivity
connecTorsion   = geo_matrix.connecTorsion; % Torsional spring connectivity
Nnodes          = geo_matrix.Nnodes;        % Total number of nodes
connecNodes     = geo_matrix.connecNodes;   % General node connectivity

%% 6. Compute Spring Constants
% Linear spring stiffness (limiting spring)
k_lim = E * (t_lim * UC) ./ (Ch_L + t + t_ch);
k_s   = 5 * k_lim;         % Rigid Spring stiffness

% Torsional spring stiffness
D  = (E * t_lim^3) / (12 * (1 - nu^2));
kt = D * (UC ./ (Unit_sep + t + t_ch));

% Combine linear & torsional springs
Ks = [k_s', k_lim', kt'];

% Bistable spring stiffness (nonlinear)
kBS = (E ./ R.^2) .* ( ...
    - (1.97 .* R .* t.^5) ./ Hv.^3 + ...
      (7.4  .* Hv.^2 .* t.^3) ./ R.^2 + ...
      (3.5  .* R .* t.^4) ./ Hv.^2 + ...
      (0.37 .* Hv.^2 .* t.^2) ./ R   + ...
      (42.2 .* t.^5) ./ Hv.^2         - ...
      (35.8 .* Hv .* t.^4) ./ R.^2    + ...
      (71.8 .* t.^5) ./ (Hv .* R)     - ...
      (3.4  .* R .* t.^3) ./ Hv       - ...
      (67.7 .* t.^4) ./ Hv            + ...
      4.2  .* R .* t.^2               + ...
      11.1 .* t.^3);

% Bistable spring shape parameter alpha
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

% Tip displacement
tip_dis = 0.2565 * t + 2.1425 * Hv;

% Ensure non-negative offsets
dd = tip_dis - Hv - Unit_sep - t_ch/2 - t/2;
dd(dd < 0) = 0.01;

% Combine bistable parameters (non-linear spring)
Bis_S = [kBS', al', dd'];

%% 7. Prepare for Optimization Loop
% Determine free (non-fixed) nodes and total DOFs
NoBcNodes = 1:Nnodes;
NoBcNodes(BC_nodes) = [];     % Remove fixed nodes
numVar     = length(NoBcNodes) * 2;  % 2 DOFs per free node

% Create main figure for plotting all cases
figure('Position', [300 300 1000 800]);
hold on;

%% 8. Loop Over Objects & Solve Equilibrium
for i = 1:length(Obj_sizeV)
    % Pack object info for plotting
    object_info = [theta, base_L, Obj_sizeV(i)];
    
    %% 8.1 Generate Initial Guess
    which_state = which_stateV(i);
    % Approximate deformed positions for inverted state
    node_locations = approx_full_inverted_state( ...
        which_state, number_of_segments, segment_length, ...
        segment_height, start_angle, 1.0 * dd');
    node_deformations   = node_locations - coords0;
    deformed_deltas     = reshape(node_deformations(:, NoBcNodes), [numVar, 1]);
    x0s = deformed_deltas;  % Initial guess for fmincon

    %% 8.2 Set fmincon Options & Objective
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'MaxFunctionEvaluations', 1e6, ...
        'OptimalityTolerance', 1e-6, ...
        'MaxIterations', 1e5, ...
        'StepTolerance', 1e-16, ...
        'SpecifyObjectiveGradient', true);

    % Objective: total elastic energy given DOFs
    Tot_En = @(sym_dof) sym_elasticenergy( ...
        sym_dof, Ks, Bis_S, coords0, connecRigid, ...
        connecBistable, connecTorsion, BC_nodes, Nnodes);

    %% 8.3 Run Optimization
    [x_sol, fval, exitflag, ~, ~, grad, hess] = ...
        fmincon(Tot_En, x0s, [], [], [], [], [], [], [], options);

    % Evaluate Hessian eigenvalues for stability check
    lam_i      = eig(hess);
    isposdef   = all(lam_i > 0);
    min_lam    = min(lam_i);
    max_lam    = max(lam_i);

    % Display solver diagnostics
    fprintf('Energy | ||Grad|| | min(lam) | max(lam) | Stable | ExitFlag\n');
    fprintf('%8.4f | %7.2e | %8.2e | %8.2e |   %d    |   %d\n', ...
        fval, norm(grad), min_lam, max_lam, isposdef, exitflag);

    %% 8.4 Plot the Equilibrium Configuration
    subplot(2, 2, i);
    plot_solution_robot( ...
        x_sol, coords0, BC_nodes, connecBistable, ...
        connecRigid, Nnodes, number_of_segments, ...
        segment_length, segment_height, object_info);
    title(sprintf('Object Size %d', i));
end

hold off;
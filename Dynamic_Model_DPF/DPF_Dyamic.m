%{
 Name: Juan Osorio
 Date: APR 2025
 Task: Spring Lattice Model to approximate Dome Phalanx Finger
%}

%% Clear workspace and set up figure defaults
clc; close all; clear;                          % Clear command window, close figures, clear variables
set(0, 'DefaultFigureColor', 'white', ...        % White background for figures
    'DefaultAxesFontSize', 22, ...               % Large font size for readability
    'DefaultAxesFontName', 'CMU Serif', ...      % Serif font for axes labels
    'DefaultLineLineWidth', 2, ...               % Thicker lines by default
    'DefaultAxesBox', 'on');                     % Box around axes
colors = {'#440154','#453781','#2d708e','#20a387','#55c667','#b8de29'}; % Color palette for plots

%% ================= Geometric Data ====================================
%----- Finger geometry parameters -----
rb = 5;                                          % Base radius of dome units
Hv = [4.0 4.5 5.0 5.0 3.0 3.0];                  % Dome heights for each segment
Ch_L = 6.0 * ones(size(Hv));                    % Characteristic length of bistable units
Unit_sep = 1.0 * ones(size(Hv));                % Separation between adjacent units

t = 0.83;                                        % Thickness of dome shell
t_lim = 1.0;                                     % Limiter thickness

% Problem-specific parameters
t_ch = 1.0;                                      % Connector thickness between units
number_of_segments = length(Hv);                 % Number of dome segments in the finger

% Material properties (MPa)
E  = 26.0;                                       % Young's modulus of material
nu = 0.35;                                       % Poisson's ratio of material

% Derived geometric quantities
base = 2 * rb;                                  % Full base diameter
UC   = 1.5 * base;                              % Unit cell width
R    = (Hv.^2 + rb^2) ./ (2 * Hv);              % Radius of curvature for each dome
start_angle = 0;                                 % Initial angular offset for geometry
BC_nodes = [1 2 3];                              % Fixed boundary nodes for constraints

% Segment dimensions for lattice generation
segment_length = Ch_L + t_ch + Unit_sep;        % Effective length of each segment
segment_height = UC/2 + t_lim + t;               % Effective height of each segment

% Generate nodal coordinates and connectivity via helper function
geo_matrix = DPG_Geometry(segment_height, segment_length, number_of_segments);
coords0       = geo_matrix.coords0;               % Reference (undeformed) node coordinates
connecRigid   = geo_matrix.connecRigid;           % Rigid connections (pinned joints)
connecBistable= geo_matrix.connecBistable;        % Bistable spring connections
connecTorsion = geo_matrix.connecTorsion;         % Torsional spring connections
Nnodes        = geo_matrix.Nnodes;                % Total number of nodes in mesh
connecNodes   = geo_matrix.connecNodes;           % Full connectivity for dynamic system


% ================Spring Data==============================================
k_lim = E*(t_lim*UC)./(Ch_L+t+t_ch); % Linear Spring
k_s = 5*k_lim;

D = (E*(t_lim)^3)/(12*(1-nu^2));
kt = D*(UC./(Unit_sep+t+t_ch));                           % Torsional Spring

Ks = [k_s',k_lim',kt'];

% Bistable Spring
kBS = (E ./ R.^2) .* ( ...
    - (1.97 .* R .* t.^5) ./ Hv.^3 + ...
    (7.4 .* Hv.^2 .* t.^3) ./ R.^2 + ...
    (3.5 .* R .* t.^4) ./ Hv.^2 + ...
    (0.37 .* Hv.^2 .* t.^2) ./ R + ...
    (42.2 .* t.^5) ./ Hv.^2 - ...
    (35.8 .* Hv .* t.^4) ./ R.^2 + ...
    (71.8 .* t.^5) ./ (Hv .* R) - ...
    (3.4 .* R .* t.^3) ./ Hv - ...
    (67.7 .* t.^4) ./ Hv + ...
    4.2 .* R .* t.^2 + ...
    11.1 .* t.^3);

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

tip_dis =  0.2565 * t + 2.1425 * Hv;

dd = tip_dis - Hv - Unit_sep - t_ch/2 - t/2;

for i = 1:length(dd)
    if dd(i)<0
        dd(i) = 0.01;
    end
end

Bis_S = [kBS',al',dd'];

% ================Dynamic Data ============================================
eta_iso = 0.01;
eta = 0.05;

% =========================== Spatial Symbolic Setup=======================
numVar = (Nnodes-length(BC_nodes))*2;

%===========================Dynamic Setup=========================
% PRESSURE PARAMETERS
P0 = 0.7; % Presure [MPa]
k_dome = 1.0;

tau_1 = 0.0;
tau_2 = 1.0;
tau_3 = 5.0;
tau_4 = 5.1;


% MASS PARAMETERS
DPG_mass = (7*10^-6)*(number_of_segments/5); %[Ton]
mass_V = (DPG_mass/Nnodes)*ones(numVar,1);

% =========================================================================
% Integration scheme for time (RK45)
% =========================================================================
%===========================Initial State=========================
which_state = 1;

NoBcNodes = 1:Nnodes;
NoBcNodes(BC_nodes) = [];
numVar = length(NoBcNodes)*2;

node_locations_deformed = approx_full_inverted_state(which_state,number_of_segments,segment_length,segment_height,start_angle,1.2*dd');
node_deformations = node_locations_deformed - coords0;
deformed_coord_deltas = reshape(node_deformations(:,NoBcNodes),[numVar,1]);

x0s = deformed_coord_deltas;

options = optimoptions('fmincon','Algorithm','sqp', ...
    'MaxFunctionEvaluations',1000000,'OptimalityTolerance',1.0000e-6,'MaxIterations',100000,...
    'StepTolerance',1.0000e-16,'SpecifyObjectiveGradient',true);

Tot_En = @(sym_dof) sym_elasticenergy(sym_dof,Ks,Bis_S,coords0,connecRigid,connecBistable,connecTorsion,BC_nodes,Nnodes);

[x_sol,fval,exitflag,~,~,grad,hess]  = fmincon(Tot_En,x0s,[],[],[],[],[],[],[],options);

Evals = fval;
lam_i = eig(hess);
isposdef = all(lam_i > 0);
min_lam = min(lam_i);
max_lam = max(lam_i);

fprintf('Energy| norm(Grad) | min(lam)| max(lam)|Stable|Exit Flag| \n')
fprintf('%10.4f|%10.4d  |%10.4d |%10.4d |%10d |%10d |\n',Evals,norm(grad),min_lam,max_lam,isposdef,exitflag)

plot_solution(x_sol,coords0,BC_nodes,connecBistable,connecRigid,Nnodes,number_of_segments,segment_length,segment_height)
%===========================RUN ANALYSIS=========================
%x0s = zeros(numVar,1);
x0s = x_sol;
v0 = zeros(numVar,1);
x0 = [x0s;v0];

TT = 15.0;
tspan = [0, TT];
[time,x] = ...
    solve_system(x0,tspan,numVar,mass_V,eta,eta_iso,k_dome,P0,tau_1,tau_2,tau_3, tau_4,Ks,Bis_S,connecNodes,connecRigid,connecBistable,connecTorsion,BC_nodes,coords0,Nnodes);

display('System Solved!')

%% =========================== Time dependent Pressure ============
P_time = zeros(length(time),1);

for i = 1:length(time)
    P_time(i) = P_t(P0,time(i), tau_1, tau_2, tau_3, tau_4);
end

%% Calculate Energy and Fin
E_t = zeros(length(time),1);
F_in_t = zeros(length(time),1);
F_d_in_t = zeros(length(time),1);
F_d_iso = zeros(length(time),1);

for i = 1:length(time)
    E_tot = sym_elasticenergy(x(i,1:numVar)',Ks,Bis_S,coords0,connecRigid,connecBistable,connecTorsion,BC_nodes,Nnodes);
    grad_E = Grad_sym_elasticenergy(x(i,1:numVar)',Ks,Bis_S,coords0,connecRigid,connecBistable,connecTorsion,BC_nodes,Nnodes);
    E_t(i) = E_tot;
end

%% =========================================================================
% Plot final state
% =========================================================================
plot_solution(x(end,1:numVar),coords0,BC_nodes,connecBistable,connecRigid,Nnodes,number_of_segments,segment_length,segment_height)

%% =========================================================================
% Plot Energy vs time
% =========================================================================
tip_dis = x(:,numVar-1:numVar);
tip_vel = x(:,end-1:end);
tip_dis_plot = sqrt(sum(tip_dis.^2,2));

% Energy amnd pressure Plot
figure
yyaxis left
plot(time,E_t,'-','Color',colors{3},'LineWidth', 3);hold on
ylabel('Strain Energy [mJ] (Solid Line)', 'FontWeight', 'bold')
ax = gca;
ax.YColor = 'k';  % Change to your desired color

yyaxis right
plot(time,P_time/P0,'--','Color',colors{1},'LineWidth', 3)
ylabel('Normalized P [-] (Dashed Line)');
ylim([0 2]);
ax = gca;
ax.YColor = 'k';  % Change to your desired color
xlabel('Time [s]');


%%
%=========================================================================
%Animation
%=========================================================================
Dynamic_Animation('Pressure_Dynamic.gif',x,coords0,BC_nodes,connecBistable,connecRigid,Nnodes,numVar,number_of_segments,segment_length,segment_height)

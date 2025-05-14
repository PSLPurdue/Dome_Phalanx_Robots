%{
Name: Juan Osorio
Date: Apr 2025
Task: Spring Lattice Model to approximate Dome Phalanx Finger
%}

clc; close all; clear;

% ================Geometric Data============================================
% FINGER PARAMETERS
rb = 5;

Hv = [5.0 5.0 3.0];
Ch_L = 6.5*ones(size(Hv));
Unit_sep = 1.0*ones(size(Hv));

t = 0.8;
t_lim = 1.25;

% PROBLEM PARAMETER
t_ch = 1.0;
number_of_segments = length(Hv);

E =  26.0;
nu = 0.35;

base = 2*rb;
UC = 1.5*base;
R = (Hv.^2 + rb^2)./(2*Hv);

start_angle = 0;
BC_nodes = [1 2 3];

segment_length = Ch_L + t_ch + Unit_sep;
segment_height = UC/2 + t_lim + t;

geo_matrix =  DPG_Geometry(segment_height,segment_length,number_of_segments);
coords0 = geo_matrix.coords0;
connecRigid = geo_matrix.connecRigid;
connecBistable = geo_matrix.connecBistable;
connecTorsion = geo_matrix.connecTorsion;
Nnodes = geo_matrix.Nnodes;
connecNodes = geo_matrix.connecNodes;

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

al = 1.3*(((0.2941 .* Hv.^3) ./ R.^3 - ...
    (0.4412 .* t.^3) ./ Hv.^3 - ...
    (2.2059 .* Hv.^2 .* t) ./ R.^3 - ...
    (0.7353 .* Hv.^2) ./ R.^2 - ...
    (6.9118 .* t.^3) ./ (Hv.^2 .* R) + ...
    (1.1765 .* t.^2) ./ Hv.^2 - ...
    (3.0882 .* t.^3) ./ (Hv .* R.^2) + ...
    (7.0588 .* Hv .* t) ./ R.^2 + ...
    (13.3824 .* t.^2) ./ (Hv .* R) + ...
    (0.4412 .* Hv) ./ R + ...
    (13.6765 .* t.^3) ./ R.^3 - ...
    (15.0 .* t.^2) ./ R.^2 - ...
    (4.7059 .* t) ./ R));

tip_dis =  0.2565 * t + 2.1425 * Hv;

dd = tip_dis - Hv - Unit_sep - t_ch/2 - t/2;

for i = 1:length(dd)
    if dd(i)<0
        dd(i) = 0.01;
    end
end

Bis_S = [kBS',al',dd'];

% ===========================Symbolic Setup=========================
NoBcNodes = 1:Nnodes;
NoBcNodes(BC_nodes) = [];
numVar = length(NoBcNodes)*2;


%% ===========================Initial Guess ================================
which_state = 2^number_of_segments;

node_locations_deformed = approx_full_inverted_state(which_state,number_of_segments,segment_length,segment_height,start_angle,dd');
node_deformations = node_locations_deformed - coords0;
deformed_coord_deltas = reshape(node_deformations(:,NoBcNodes),[numVar,1]);

x0s = deformed_coord_deltas;

%% =========Parameters for optimization ====================================
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

%% Plot Solution
plot_solution(x_sol,coords0,BC_nodes,connecBistable,connecRigid,Nnodes,number_of_segments,segment_length,segment_height)

%plot_solution_robot(x_sol,coords0,BC_nodes,connecBistable,connecRigid,Nnodes,number_of_segments,segment_length,segment_height)

% % Set font properties
% set(gca, 'FontName', 'CMU Serif', 'FontSize', 14); % for axes
% set(findall(gcf, '-property', 'FontName'), 'FontName', 'CMU Serif'); % for all text
% set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14); % set specific font size for all text
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [0 0 3.5 3.5]); % [left bottom width height]
% % Save the figure
% print(strcat('Energy_Solution'), '-dsvg', '-r300'); % save as PNG at 300 DPI

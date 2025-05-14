% run_DPW_model_variant.m
% Name: Juan Osorio
% Date: Apr 2025
% Task: Spring Lattice Model to approximate Dome Phalanx Finger (alternate Hs order)
%   - Simulates dynamic response of a multi-section bistable metasheet "finger"
%   - Plots normalized tip trajectory for both metastable and monostable cases
%   - Creates animated GIF of finger motion

clc; close all; clear;

%% 1. Configure figure defaults for consistent styling
set(0, 'DefaultFigureColor', 'white', ...      % white background
       'DefaultAxesFontSize', 22, ...          % font size
       'DefaultAxesFontName','CMU Serif', ...  % font family
       'DefaultLineLineWidth', 2, ...          % line width
       'DefaultAxesBox', 'on');                % boxed axes
colors = {'#440154','#453781','#2d708e','#20a387','#55c667','#b8de29'};  % color palette

%% 2. Define simulation scenarios
Hs = [1.1, 2.5];     % dome heights for phase 1
P0 = 0.4;            % peak pressure magnitude


figure; hold on;     % prepare figure for tip path plots

%% 3. Loop over each height scenario
for p_i = 1:length(Hs)
    % 3.1 Global parameters common to all sections
    t      = 1.0;      % dome thickness
    t_lim  = 1.0;      % limiting layer thickness
    TT     = 5;        % total simulation time
    tspan  = linspace(0, TT, 24000);  % time vector
    rb     = 5;        % dome base radius

    %% 3.2 Section 1: Two-unit finger base
    Hv       = [Hs(p_i), Hs(p_i)];             % phase 1 dome height
    Ch_L     = 6.5 * ones(size(Hv));          % channel length per unit
    Unit_sep = 1.0 * ones(size(Hv));          % unit-cell separation
    which_state = 1;                          % initial shape index
    tau_1 = 0; tau_2 = 0.1; tau_3 = 1.1; tau_4 = 1.2;  % pressure timing
    % Run setup, optimization, and dynamic simulation
    [time_1, x_1, connecNodes_1, Nnodes_1, BC_nodes_1, coords0_1] = ...
        section_i_robot(Hv, Ch_L, Unit_sep, t, t_lim, rb, which_state, ...
                       P0, tau_1, tau_2, tau_3, tau_4, tspan);
    
    %% 3.3 Section 2: Single-unit follower
    Hv       = [2.0];                           % follower dome height
    Ch_L     = 6.5 * ones(size(Hv));
    Unit_sep = 1.0 * ones(size(Hv));
    which_state = 1;
    tau_1 = 1.1; tau_2 = 1.2; tau_3 = 2.2; tau_4 = 2.3;  % follower pressure timing
    [time_2, x_2, connecNodes_2, Nnodes_2, BC_nodes_2, coords0_2] = ...
        section_i_robot(Hv, Ch_L, Unit_sep, t, t_lim, rb, which_state, ...
                       P0, tau_1, tau_2, tau_3, tau_4, tspan);

    %% 3.4 Section 3: 3D end-effector section
    Hv       = [Hs(p_i) + 0.5, 4.5, 4.5, 4.5];    % end section dome heights
    Ch_L     = 6.5 * ones(size(Hv));
    Unit_sep = 1.0 * ones(size(Hv));
    which_state = 1;
    [time_3, x_3, connecNodes_3, Nnodes_3, BC_nodes_3, coords0_3] = ...
        section_i_robot(Hv, Ch_L, Unit_sep, t, t_lim, rb, which_state, ...
                       P0, tau_1, tau_2, tau_3, tau_4, tspan);
    
    %% 4. Compute and store tip trajectory
    tip_indices = linspace(1, size(x_1,1)-1, 500);  % sample frames
    tip_V = zeros(length(tip_indices), 3);
    for k = 1:length(tip_indices)
        idx = floor(tip_indices(k));
        % 4.1 Reconstruct section1 coords and orientation
        free1 = setdiff(1:Nnodes_1, BC_nodes_1);
        uvw1 = zeros(2, Nnodes_1);
        uvw1(:, free1) = reshape(x_1(idx,1:2*length(free1)), [2, length(free1)]);
        coords1 = coords0_1 + uvw1;
        coords1(2,:) = -coords1(2,:);     % mirror y-axis downwards
        vec1 = coords1(:,end-1) - coords1(:,end-2);
        theta1 = -atan2(vec1(2), vec1(1));
        % 4.2 Section2 attach and orient
        free2 = setdiff(1:Nnodes_2, BC_nodes_2);
        uvw2 = zeros(2, Nnodes_2);
        uvw2(:, free2) = reshape(x_2(idx,1:2*length(free2)), [2,length(free2)]);
        coords2 = coords0_2 + uvw2;
        R2 = [cos(theta1), -sin(theta1); sin(theta1), cos(theta1)];
        coords2 = R2 * coords2;
        coords2(1,:) = coords2(1,:) + coords1(1,end-1);
        coords2(2,:) = coords2(2,:) - coords1(2,end-1);
        vec2 = coords2(:,end) - coords2(:,end-3);
        theta2 = atan2(vec2(2), vec2(1));
        % 4.3 Section3 3D rotation & translation
        free3 = setdiff(1:Nnodes_3, BC_nodes_3);
        uvw3 = zeros(2, Nnodes_3);
        uvw3(:, free3) = reshape(x_3(idx,1:2*length(free3)), [2, length(free3)]);
        coords3 = coords0_3 + uvw3;
        Rx = [1,0,0; 0,cosd(90),-sind(90); 0,sind(90),cosd(90)];
        Rz = [cos(theta1+theta2), -sin(theta1+theta2), 0; sin(theta1+theta2), cos(theta1+theta2), 0; 0,0,1];
        coords3_3D = [coords3; zeros(1,size(coords3,2))];
        coords3_3D = Rz * Rx * coords3_3D;
        coords3_3D(3,:) = -coords3_3D(3,:);   % mirror z-axis downwards
        coords3_3D = coords3_3D + [coords2(:,end); 0];
        tip_V(k,:) = coords3_3D(:,end)';
    end
    % 4.4 Normalize and plot tip path
    tip_norm = tip_V - tip_V(1,:);
    
    plot3(tip_norm(:,1), tip_norm(:,2), tip_norm(:,3), '-', ...
          'LineWidth', 3.0, 'Color', colors{p_i+1});
end

%% 5. Annotate plot with legend and axes
plot3(0, 0, 0, 'o', 'MarkerSize', 15, 'Color', colors{4}, 'MarkerFaceColor', colors{4});
legend('Monostable','Metastable', 'Start');
box on;
xlabel('Tip U [mm]'); ylabel('Tip V [mm]'); zlabel('W [mm]');
view(90,0);
hold off;

%% 6. Create animated GIF for last simulation scenario
figure;
Dynamic_Animation('DPW_H.gif', x_1, Nnodes_1, BC_nodes_1, connecNodes_1, coords0_1, ...
                  x_2, Nnodes_2, BC_nodes_2, connecNodes_2, coords0_2, ...
                  x_3, Nnodes_3, BC_nodes_3, connecNodes_3, coords0_3);

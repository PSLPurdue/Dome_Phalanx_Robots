function [] = Dynamic_Animation(name, x_1, Nnodes_1, BC_nodes_1, connecNodes_1, coords0_1, x_2, Nnodes_2, BC_nodes_2, connecNodes_2, coords0_2, x_3, Nnodes_3, BC_nodes_3, connecNodes_3, coords0_3)
% DYNAMIC_ANIMATION  Create animated GIF of multi-segment robot motion
%   Inputs:
%     name           – output GIF filename
%     x_1            – time-series solution vector for segment 1
%     Nnodes_1       – number of nodes in segment 1
%     BC_nodes_1     – fixed node indices for segment 1
%     connecNodes_1  – connectivity pairs for segment 1
%     coords0_1      – original coordinates for segment 1
%     x_2,...        – analogous inputs for segment 2
%     x_3,...        – analogous inputs for segment 3

%% 1. Figure and styling setup
set(0, 'DefaultFigureColor', 'white', ...        % white background
       'DefaultAxesFontSize', 22, ...            % larger fonts
       'DefaultAxesFontName','CMU Serif', ...    % serif font
       'DefaultLineLineWidth', 3, ...           % thick lines
       'DefaultAxesBox', 'on');                 % boxed axes
colors = {'#440154','#453781','#2d708e','#20a387','#55c667','#b8de29'}; % color palette

%% 2. Initialize animation variables
tip_V = [];                                   % store tip trajectory
anim_V = linspace(1, size(x_1,1)-1, 500);     % frames sampled uniformly
h_fig = figure;                                % create new figure

%% 3. Main animation loop over frames
for anim_i = 1:length(anim_V)
    frame_idx = floor(anim_V(anim_i));         % integer frame index
    
    %% 3.1 Section 1: Segment 1 deformation
    % 3.1.1 Identify free DOFs for segment 1
    NoBcNodes_1 = 1:Nnodes_1;
    NoBcNodes_1(BC_nodes_1) = [];
    numVar_1 = (Nnodes_1 - length(BC_nodes_1)) * 2;
    % 3.1.2 Extract solution and reshape into displacements
    dof = x_1(frame_idx, 1:numVar_1);
    uvw = zeros(2, Nnodes_1);
    uvw(:, NoBcNodes_1) = reshape(dof, [2, length(NoBcNodes_1)]);
    coords_1 = coords0_1 + uvw;               % deformed coords
    % 3.1.3 Mirror vertically (flip y-axis)
    coords_1_R = coords_1;
    coords_1_R(2,:) = -coords_1(2,:);
    % 3.1.4 Compute rotation angle based on last link
    x_ij_1 = coords_1(:, end-1) - coords_1(:, end-2);
    theta_1 = -atan2(x_ij_1(2), x_ij_1(1));     % negative for desired orientation
    % 3.1.5 Plot connections for segment 1
    for pair = connecNodes_1'
        node1 = pair(1); node2 = pair(2);
        line([coords_1_R(1,node1), coords_1_R(1,node2)], ...
             [coords_1_R(2,node1), coords_1_R(2,node2)], [0,0], ...
             'Color', colors{3}); hold on;
    end
    % 3.1.6 Plot nodes for segment 1
    for i = 1:Nnodes_1
        plot(coords_1_R(1,i), coords_1_R(2,i), 'o', ...
             'Color','k','MarkerFaceColor',colors{6}, 'MarkerSize',8);
    end
    
    %% 3.2 Section 2: Segment 2 deformation and placement
    % 3.2.1 Identify free DOFs
    NoBcNodes_2 = 1:Nnodes_2;
    NoBcNodes_2(BC_nodes_2) = [];
    numVar_2 = (Nnodes_2 - length(BC_nodes_2)) * 2;
    % 3.2.2 Extract and reshape displacements
    dof = x_2(frame_idx, 1:numVar_2);
    uvw = zeros(2, Nnodes_2);
    uvw(:, NoBcNodes_2) = reshape(dof, [2, length(NoBcNodes_2)]);
    coords_2 = coords0_2 + uvw;
    % 3.2.3 Build rotation matrix from theta_1
    R_2 = [cos(theta_1), -sin(theta_1); sin(theta_1), cos(theta_1)];
    % 3.2.4 Rotate then translate segment 2 to attach to segment 1 tip
    coords_2_T = R_2 * coords_2;
    coords_2_T(1,:) = coords_2_T(1,:) + coords_1(1,end-1);
    coords_2_T(2,:) = coords_2_T(2,:) - coords_1(2,end-1);
    % 3.2.5 Compute local rotation for segment 2
    x_ij_2 = coords_2(:,end) - coords_2(:,end-3);
    theta_2 = atan2(x_ij_2(2), x_ij_2(1));
    % 3.2.6 Plot segment 2 connections and nodes
    for pair = connecNodes_2'
        node1 = pair(1); node2 = pair(2);
        line([coords_2_T(1,node1), coords_2_T(1,node2)], ...
             [coords_2_T(2,node1), coords_2_T(2,node2)], [0,0], ...
             'Color', colors{1});
    end
    for i = 1:Nnodes_2
        plot(coords_2_T(1,i), coords_2_T(2,i), 'o', ...
             'Color','k','MarkerFaceColor',colors{6}, 'MarkerSize',8);
    end
    
    %% 3.3 Section 3: Segment 3 deformation in 3D
    % 3.3.1 Identify free DOFs
    NoBcNodes_3 = 1:Nnodes_3;
    NoBcNodes_3(BC_nodes_3) = [];
    numVar_3 = (Nnodes_3 - length(BC_nodes_3)) * 2;
    % 3.3.2 Extract and reshape displacements
    dof = x_3(frame_idx, 1:numVar_3);
    uvw = zeros(2, Nnodes_3);
    uvw(:, NoBcNodes_3) = reshape(dof, [2, length(NoBcNodes_3)]);
    coords_3 = coords0_3 + uvw;
    % 3.3.3 Define 3D rotations: 90° about x, then theta_1+theta_2 about z
    theta_x_3 = deg2rad(90);
    theta_z_3 = theta_1 + theta_2;
    Rx = [1, 0, 0; 0, cos(theta_x_3), -sin(theta_x_3); 0, sin(theta_x_3), cos(theta_x_3)];
    Rz = [cos(theta_z_3), -sin(theta_z_3), 0; sin(theta_z_3), cos(theta_z_3), 0; 0,0,1];
    R = Rz * Rx;
    % 3.3.4 Lift 2D coords to 3D plane, apply R, then mirror z-axis
    coords_3_3D = [coords_3; zeros(1, size(coords_3,2))];
    coords_3_T = R * coords_3_3D;
    coords_3_T(3,:) = -coords_3_T(3,:);
    % 3.3.5 Translate to attach at segment 2 tip
    coords_3_T = coords_3_T + [coords_2_T(:,end); 0];
    % 3.3.6 Plot 3D connections and nodes
    for pair = connecNodes_3'
        node1 = pair(1); node2 = pair(2);
        line([coords_3_T(1,node1), coords_3_T(1,node2)], ...
             [coords_3_T(2,node1), coords_3_T(2,node2)], ...
             [coords_3_T(3,node1), coords_3_T(3,node2)], 'Color', colors{1});
    end
    for i = 1:Nnodes_3
        plot3(coords_3_T(1,i), coords_3_T(2,i), coords_3_T(3,i), 'o', ...
             'Color','k','MarkerFaceColor',colors{6}, 'MarkerSize',8);
    end
    % 3.3.7 Accumulate and plot tip trajectory
    tip_V = [tip_V; coords_3_T(:,end)'];
    plot3(tip_V(:,1), tip_V(:,2), tip_V(:,3), '-o', 'LineWidth',1.5, 'MarkerSize',4);

    %% 3.4 Finalize frame and save to GIF
    hold off;
    view([1 1 1]);            % isometric view
    xlabel('X'); ylabel('Y'); zlabel('Z');
    xlim([0 70]); ylim([-30 30]); zlim([-30 20]);
    make_animation(h_fig, anim_i, name);
    clf;                      % clear figure for next frame
    pause(0.01);             % brief pause for rendering
end


%% 4. Helper: Write or append frame to GIF file
    function make_animation(h, index, filename)
        drawnow;
        frame = getframe(h);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if index == 1
            % first frame: create new GIF
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
        else
            % subsequent frames: append
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
        end
    end
end

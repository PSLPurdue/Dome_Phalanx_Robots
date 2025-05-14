function [] = plot_solution_robot(state_Vec, coords0, BC_nodes, connecBistable, connecRigid, Nnodes, number_of_segments, segment_length, segment_height, object_info)
% PLOT_SOLUTION_ROBOT  Visualize equilibrium configuration and its reflection
%   Given the optimized state vector, original node coords, and spring
%   connectivity, plot the deformed spring lattice and its mirror image
%   about a horizontal line at y = base_L.
%
%   Inputs:
%     state_Vec         – [1×2*(Nnodes-#BC)] optimized displacements for free nodes
%     coords0           – [2×Nnodes] undeformed node coordinates
%     BC_nodes          – indices of fixed (boundary) nodes
%     connecBistable    – [Q×2] bistable spring node pairs
%     connecRigid       – [2P×2] rigid & limit-layer spring node pairs
%     Nnodes            – total number of nodes
%     number_of_segments– number of segments in the “finger” chain
%     segment_length    – [1×Nsegments] original segment lengths
%     segment_height    – height of each rectangular segment
%     object_info       – [theta, base_L, Obj_size] for rotation & base lines

%% 1. Unpack object_info and reset figure defaults
theta    = object_info(1);  % global rotation angle [deg]
base_L   = object_info(2);  % y-coordinate of mirror axis
Obj_size = object_info(3);  % half-height of object to visualize

% Ensure consistent styling for plots
set(0, 'DefaultFigureColor', 'white', ...
    'DefaultAxesFontSize', 22, ...
    'DefaultAxesFontName','CMU Serif', ...
    'DefaultLineLineWidth', 2, ...
    'DefaultAxesBox', 'on');
colors = {'#440154','#453781','#2d708e','#20a387','#55c667','#b8de29'};
hold on;

%% 2. Compute deformed node coordinates
% 2.1 Map free DOFs back into full UVW displacement array
NoBcNodes = 1:Nnodes;
NoBcNodes(BC_nodes) = [];
uvw = zeros(2, Nnodes);
uvw(:, NoBcNodes) = reshape(state_Vec, [2, length(NoBcNodes)]);
coords = coords0 + uvw;  % apply displacements

% 2.2 Shift downward by segment_height for correct baseline
coords(2,:) = coords(2,:) - segment_height;

% 2.3 Rotate all nodes by theta about the origin
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
for i = 1:Nnodes
    coords(:,i) = R * coords(:,i);
end

%% 3. Plot original (bottom) lattice
% 3.1 Bistable springs: color by extension (active vs. unactive)
for pair = connecBistable'
    i = pair(1); j = pair(2);
    % length change comparison
    rest_len = norm(coords0(:,i)-coords0(:,j));
    cur_len  = norm(coords(:, i) - coords(:, j));
    if round(cur_len - rest_len) > 0
        col = colors{3};  % extended
    else
        col = colors{1};  % compressed or neutral
    end
    % draw line
    line([coords(1,i) coords(1,j)], [coords(2,i) coords(2,j)], ...
         'Color', col, 'LineWidth', 2);
end

% 3.2 Rigid & limit-layer springs: always black
for pair = connecRigid'
    i = pair(1); j = pair(2);
    line([coords(1,i) coords(1,j)], [coords(2,i) coords(2,j)], ...
         'Color', 'k', 'LineWidth', 2);
end

% 3.3 Nodes as filled circles
for i = 1:Nnodes
    plot(coords(1,i), coords(2,i), 'o', 'MarkerSize', 10, ...
         'Color', 'k', 'MarkerFaceColor', colors{6});
end

% 3.4 Draw base and object-size reference lines
yline(base_L,       '--k', 'LineWidth', 2);
yline(base_L - Obj_size/2, '--', 'Color', colors{1}, 'LineWidth', 2);
yline(base_L + Obj_size/2, '--', 'Color', colors{1}, 'LineWidth', 2);

%% 4. Plot mirrored (top) lattice about y = base_L
% 4.1 Reflect coords vertically
coords_mir = coords;
coords_mir(2,:) = 2*base_L - coords_mir(2,:);

% 4.2 Bistable springs (same color logic)
for pair = connecBistable'
    i = pair(1); j = pair(2);
    rest_len = norm(coords0(:,i)-coords0(:,j));
    cur_len  = norm(coords_mir(:, i) - coords_mir(:, j));
    if round(cur_len - rest_len) > 0
        col = colors{3};
    else
        col = colors{1};
    end
    line([coords_mir(1,i) coords_mir(1,j)], [coords_mir(2,i) coords_mir(2,j)], ...
         'Color', col, 'LineWidth', 2);
end

% 4.3 Rigid springs mirrored
for pair = connecRigid'
    i = pair(1); j = pair(2);
    line([coords_mir(1,i) coords_mir(1,j)], [coords_mir(2,i) coords_mir(2,j)], ...
         'Color', 'k', 'LineWidth', 2);
end

% 4.4 Plot mirrored nodes
for i = 1:Nnodes
    plot(coords_mir(1,i), coords_mir(2,i), 'o', 'MarkerSize', 10, ...
         'Color', 'k', 'MarkerFaceColor', colors{6});
end

%% 5. Final plot adjustments
xlabel('x [mm]');  ylabel('y [mm]');
% Expand axis to show both halves
xlim([-10, (number_of_segments+1)*mean(segment_length)+3]);
ylim([min(min(coords(2,:)),min(coords_mir(2,:)))-5, ...
      max(max(coords(2,:)),max(coords_mir(2,:)))+5]);

% Legend illustrating active vs. inactive bistable units
h(1) = plot(NaN,NaN,'-','Color',colors{1});
h(2) = plot(NaN,NaN,'-','Color',colors{3});
legend(h, {'Unactive Unit','Active Unit'}, 'Location','east');

hold off;
end
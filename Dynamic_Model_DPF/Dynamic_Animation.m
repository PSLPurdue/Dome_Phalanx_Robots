function [] = Dynamic_Animation(name, state_Vec, coords0, BC_nodes, connecBistable, connecRigid, Nnodes, numVar, number_of_segments, segment_length, segment_height)
% DYNAMIC_ANIMATION  Animate 2D deformation over time for a single section
%   Generates a GIF showing bistable and rigid springs deforming according
%   to `state_Vec`, highlighting active bistable springs by color.
%
%   Inputs:
%     name             – filename for output GIF
%     state_Vec        – [T×numVar] time-series displacement states
%     coords0          – [2×Nnodes] original node coordinates
%     BC_nodes         – indices of fixed (boundary) nodes
%     connecBistable   – [Q×2] bistable spring connectivities
%     connecRigid      – [P×2] rigid spring connectivities
%     Nnodes           – total number of nodes
%     numVar           – number of displacement DOFs (2*(Nnodes-#BC))
%     number_of_segments – number of segments in this section
%     segment_length   – [1×N] original lengths of each segment
%     segment_height   – height of each segment (uniform)

%% 1. Configure figure styling
set(0, 'DefaultFigureColor', 'white', ...      % white background
       'DefaultAxesFontSize', 22, ...          % axes font size
       'DefaultAxesFontName','CMU Serif', ...  % serif font
       'DefaultLineLineWidth', 3, ...          % thicker lines
       'DefaultAxesBox', 'on');                % draw box around axes
colors = {'#440154','#453781','#2d708e','#20a387','#55c667','#b8de29'};  % palette

%% 2. Determine free-node indices and animation frames
NoBcNodes = 1:Nnodes;
NoBcNodes(BC_nodes) = [];                     % remove fixed nodes
anim_V = linspace(1, size(state_Vec,1)-1, 100); % sample 100 frames uniformly
h_fig = figure;                              % create figure handle

%% 3. Loop through sampled frames and plot
for idx = 1:length(anim_V)
    frame = floor(anim_V(idx));              % current frame index

    % 3.1 Extract and reshape DOFs into nodal displacements
    dof = state_Vec(frame, 1:numVar);
    uvw = zeros(2, Nnodes);
    uvw(:, NoBcNodes) = reshape(dof, [2, length(NoBcNodes)]);
    coords = coords0 + uvw;                  % deformed coordinates

    hold on;
    % 3.2 Plot bistable springs, colored by extension/compression
    for pair = connecBistable'
        i = pair(1); j = pair(2);
        restVec = coords0(:, i) - coords0(:, j);
        restLen = norm(restVec);
        curVec  = coords(:, i) - coords(:, j);
        curLen  = norm(curVec);
        % choose color: extended vs. compressed
        if round(curLen - restLen) > 0
            col = colors{3};  % extended color
        else
            col = colors{1};  % compressed color
        end
        % draw spring
        line([coords(1,i), coords(1,j)], [coords(2,i), coords(2,j)], ...
             'Color', col);
    end

    % 3.3 Plot rigid springs in black
    for pair = connecRigid'
        i = pair(1); j = pair(2);
        line([coords(1,i), coords(1,j)], [coords(2,i), coords(2,j)], ...
             'Color', 'k');
    end

    % 3.4 Plot nodes as filled circles
    for n = 1:Nnodes
        plot(coords(1,n), coords(2,n), 'o', 'MarkerSize',10, ...
             'Color','k','MarkerFaceColor',colors{6});
    end

    % 3.5 Finalize plot limits and labels
    xlabel('x [mm]'); ylabel('y [mm]');
    xlim([0, (number_of_segments+1) * mean(segment_length) + 3]);
    ylim([-10, number_of_segments * segment_height]);
    hold off;

    % 3.6 Capture and append frame to GIF
    make_animation(h_fig, idx, name);
    clf; pause(0.01);
end

%% 4. Nested helper: write frames to GIF
    function make_animation(h, frameIdx, filename)
        drawnow;
        frame = getframe(h);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if frameIdx == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
        end
    end
end

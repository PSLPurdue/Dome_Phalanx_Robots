function [] = plot_initial(state_Vec,coords0,BC_nodes,connecBistable,connecRigid,Nnodes,number_of_segments,segment_length,segment_height)

set(0, 'DefaultFigureColor', 'white', 'DefaultAxesFontSize', 22 ...
    , 'DefaultAxesFontName','CMU Serif', 'DefaultLineLineWidth', 3, ...
    'DefaultAxesBox', 'on');
colors = {'#440154','#453781','#2d708e','#20a387','#55c667','#b8de29'};


figure
hold on;

% Calculate new coordinates from solution vector:
NoBcNodes = 1:Nnodes;
NoBcNodes(BC_nodes) = [];

dof = state_Vec;
uvw = zeros(2,Nnodes);
uvw(:,NoBcNodes) = reshape(dof,[2,length(NoBcNodes)]);
coords = coords0 + uvw;

for pair = connecBistable'
    node1 = pair(1); node2 = pair(2);
    x = [coords(1, node1), coords(1, node2)];
    y = [coords(2, node1), coords(2, node2)];

    XD0 = coords0(:,node1) - coords0(:,node2);
    lijB0 = sqrt(sum(XD0.^2))';
    XD = coords(:,node1) - coords(:,node2);
    lijB = sqrt(sum(XD.^2))';

    if (round(lijB - lijB0) > 0)
        bi_spring_color = colors{3};
    else
        bi_spring_color = colors{1};
    end

    line(x,y, 'Color','k')
end

for pair = connecRigid'
    node1 = pair(1); node2 = pair(2);
    x = [coords(1, node1), coords(1, node2)];
    y = [coords(2, node1), coords(2, node2)];
    line(x,y,'Color','black')
end

for i = 1:Nnodes
    plot(coords(1,i),coords(2,i),'o','Color','k','MarkerFaceColor',colors{6}, 'MarkerSize',10);
end

hold off;
xlabel('x [mm]'); ylabel('y [mm]');
xlim([0, (number_of_segments+1)*mean(segment_length) + 3]);
ylim([-10, (number_of_segments)*segment_height]);

end
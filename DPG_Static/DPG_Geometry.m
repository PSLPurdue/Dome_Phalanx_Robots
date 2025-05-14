% Function: DPG_Geometry
% Builds the geometric connectivity for a dome-phalanx finger (DPG) model
%
% Syntax:
%   geo_matrix = DPG_Geometry(segment_height, segment_length, number_of_segments)
%
% Inputs:
%   segment_height      - Height of each segment
%   segment_length      - [1×N] vector of lengths for each segment
%   number_of_segments  - Total number of segments
%
% Output:
%   geo_matrix - Struct containing:
%     coords0        - Initial node coordinates
%     connecRigid    - Pairs of nodes for rigid springs
%     connecBistable - Pairs of nodes for bistable springs
%     connecTorsion  - Triplets of nodes for torsional springs
%     Nnodes         - Total number of nodes
%     connecNodes    - Full connectivity list

function geo_matrix =  DPG_Geometry(segment_height,segment_length,number_of_segments)


% ============= Conectivity Matrix ========================================
connecNodes = [1 2]; %connectivity matrix
coords0 = [0    0;
    0    segment_height] ;%coordinates ( row:x,y ; columns: node#)
connecTorsion = [];

segment_length_i = 0;
for i = 1:number_of_segments

    %populate connectivity matrix
    nextnodes = max(connecNodes) + [1, 2];
    for n = nextnodes
        connecNodes = [connecNodes; n n-2];
    end
    connecNodes = [connecNodes; nextnodes];

    %matrix of points used to calculate thetas
    %      points = [1 3 2; 2 1 4; 3 4 1; 4 2 3] + 2*(i-1); %vertex point, leg point, leg point
    points = [2 1 4; 4 2 3] + 2*(i-1); %vertex point, leg point, leg point
    connecTorsion = [connecTorsion; points];

    segment_length_i = segment_length_i + segment_length(i);

    %populate coordinate matrix
    nextcoords = [0    0; 0    segment_height] + [segment_length_i; 0];
    coords0 = [coords0, nextcoords];
end


[Nsprings, ~] = size(connecNodes);
connecRigid = [];
connecBistable = [];

for i = 1:Nsprings
    connection = connecNodes(i,:);
    if (rem(connection(1),2)~=0 && rem(connection(2),2)~=0)
        connecBistable = [connecBistable; connection];
    else
        connecRigid = [connecRigid; connection];
    end
end

temp_connecBistable = connecBistable;
temp_connecBistable((connecBistable==3)) = 4;
temp_connecBistable((connecBistable==4)) = 3;
connecBistable = temp_connecBistable;

temp_connecRigid = connecRigid;
temp_connecRigid((connecRigid==3)) = 4;
temp_connecRigid((connecRigid==4)) = 3;
connecRigid = temp_connecRigid;

temp_connecNodes = connecNodes;
temp_connecNodes((connecNodes==3)) = 4;
temp_connecNodes((connecNodes==4)) = 3;
connecNodes = temp_connecNodes;

temp_connecTorsion = connecTorsion;
temp_connecTorsion((connecTorsion==3)) = 4;
temp_connecTorsion((connecTorsion==4)) = 3;
connecTorsion = temp_connecTorsion;

temp_coords0 = coords0;
temp_coords0(:,3) = coords0(:,4);
temp_coords0(:,4) = coords0(:,3);
coords0 = temp_coords0;

Nnodes = max(connecNodes,[],'all');         %number of nodes


% Fun Outputs =============================================================
connecTorsion(1,:) = [];

geo_matrix.coords0 = coords0;
geo_matrix.connecRigid = connecRigid;
geo_matrix.connecBistable = connecBistable;
geo_matrix.connecTorsion = connecTorsion;
geo_matrix.Nnodes = Nnodes;
geo_matrix.connecNodes = connecNodes;

end
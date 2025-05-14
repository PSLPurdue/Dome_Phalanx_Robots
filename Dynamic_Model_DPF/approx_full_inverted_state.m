function node_locations_deformed = approx_full_inverted_state(which_state,number_of_segments,segment_length,segment_height,start_angle,extension_vector)

% APPROX_FULL_INVERTED_STATE  Approximate deformed node positions
%   Given a binary “inversion” state index, build a piecewise‐linear
%   chain of segments (“finger”) with local rotations and extensions,
%   and return the flattened [x; y] coordinates for all nodes.
%
%   Inputs:
%     which_state        – scalar index into the 2^N truth–table
%     number_of_segments – N, the number of segments in the finger
%     segment_length     – [1×N] vector of unextended segment lengths
%     segment_height     – scalar height of each segment rectangle
%     start_angle        – initial global rotation angle [deg]
%     extension_vector   – [1×N] extension amount for each segment
%
%   Output:
%     node_locations_deformed – [2×(4*N)] matrix of node coords:
%                               [x1 x2 x3 x4 | x5 x6 …] etc.

%% 1. Build binary “truth–table” of all 2^N activation patterns
bits           = number_of_segments;
numConditions  = 2^bits;            % total rows in truth table
truthtable_vector = zeros(numConditions, bits);

for col = 1:bits
    % build repeating 0/1 column pattern of length numConditions
    block = [zeros(2^(col-1),1); ones(2^(col-1),1)];
    % tile the block to fill full column
    tiled = repmat(block, numConditions / (2^col), 1);
    % assign to truth table, reverse order so bit1 is least significant
    truthtable_vector(:, bits-col+1) = tiled;
end

% extract the activation (0 = no extension, 1 = apply extension)
activation_sequence = truthtable_vector(which_state, :)';

%% 2. Define helper functions for rotation & segment coords
% 2×2 rotation matrix for angle phi (degrees)
rot_matrix_2x2 = @(phi) [cosd(phi) -sind(phi);
                        sind(phi)  cosd(phi)];

% Given origin [x0,y0], angle, and lengths, return 4 corner points
xy = @(origin, angle, len) ...
    origin + [0, 0;
              len, 0;
              len, -segment_height;
              0, -segment_height] ...
             * rot_matrix_2x2(angle);
         
%% 3. Initialize arrays for building the deformed chain
origins(1,:) = [0, 0];               % start at global origin
angles(1)    = start_angle;          % initial orientation
length_changes = activation_sequence .* extension_vector;  % extension/no

%% 4. Loop over each segment to build node dictionary & chain
for seg = 1:number_of_segments
    % 4.1 Local corner coordinates before extension
    P = xy(origins(seg,:), angles(seg), segment_length(seg));
    node_dictionary(:,:,seg) = P;

    % 4.2 Apply extension to the third corner
    extra = [length_changes(seg), 0] * rot_matrix_2x2(angles(seg));
    node_dictionary(3,:,seg) = node_dictionary(3,:,seg) + extra;

    % 4.3 Enforce connectivity: link bottom-left to previous top-right
    if seg ~= 1
        node_dictionary(4,:,seg) = node_dictionary(3,:,seg-1);
    end

    % 4.4 Update next origin & compute local rotation
    origins(seg+1,:) = node_dictionary(2,:,seg);
    local_rotations(seg) = -atand(length_changes(seg)/segment_height);
    angles(seg+1) = angles(seg) + local_rotations(seg);

    % 4.5 Subtract global origin of first segment, then flatten
    %    This aligns all segments relative to the first corner
    n = node_dictionary(:,:,seg) - node_dictionary(4,:,1);
    coord(:, 2*seg-1:2*seg+2) = [ ...
        n(4,1), n(1,1), n(3,1), n(2,1); ... % x-coords
        n(4,2), n(1,2), n(3,2), n(2,2)   ... % y-coords
    ];
end

%% 5. Fix ordering quirk: swap columns 3 & 4 to match expected layout
temp    = coord;
temp(:,3) = coord(:,4);
temp(:,4) = coord(:,3);
coord   = temp;

%% 6. Output the deformed node locations
node_locations_deformed = coord;

end
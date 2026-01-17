function nodal_points = find_nodal_points_geometric(arcs_A, arcs_B, H_R, v_axis, v_vert)
% FIND_NODAL_POINTS_GEOMETRIC - Find nodal points using vault symmetry
%
% Theory (from homework):
%   - Diagonal arcs are NOT conics (they are parallel space curves)
%   - Arcs a_i and b_j are symmetric wrt vertical plane π_ij perpendicular to axis
%   - Nodal point N_ij lies on this symmetry plane
%   - In rectified space: symmetry planes become vertical lines
%
% Method:
%   1. Transform arc points to rectified space using H_R
%   2. Find where arcs actually cross (segment intersection)
%   3. The crossing point lies on the symmetry plane by construction
%   4. Transform back to original image coordinates

nodal_points = [];

if isempty(arcs_A) || isempty(arcs_B)
    fprintf('No arcs provided.\n');
    return;
end

if nargin < 3 || isempty(H_R)
    H_R = eye(3);
end
H_R_inv = inv(H_R);

fprintf('Finding nodal points geometrically (segment intersection in rectified space)...\n');

%% Transform arcs to rectified space
arcs_A_rect = cell(size(arcs_A));
arcs_B_rect = cell(size(arcs_B));

for i = 1:length(arcs_A)
    pts = arcs_A{i};
    pts_hom = [pts, ones(size(pts,1), 1)]';
    pts_rect = H_R * pts_hom;
    pts_rect = pts_rect(1:2, :) ./ pts_rect(3, :);
    arcs_A_rect{i} = pts_rect';
end

for i = 1:length(arcs_B)
    pts = arcs_B{i};
    pts_hom = [pts, ones(size(pts,1), 1)]';
    pts_rect = H_R * pts_hom;
    pts_rect = pts_rect(1:2, :) ./ pts_rect(3, :);
    arcs_B_rect{i} = pts_rect';
end

%% Find intersections using polyline segment crossing
% In rectified space, if arcs cross, they do so at nodal points
% The crossing defines the symmetry plane

for i = 1:length(arcs_A_rect)
    pts_A = arcs_A_rect{i};

    for j = 1:length(arcs_B_rect)
        pts_B = arcs_B_rect{j};

        % Find segment intersection
        [pt_rect, found] = find_polyline_intersection(pts_A, pts_B);

        if found
            % Transform back to original image coordinates
            pt_orig = H_R_inv * [pt_rect(:); 1];
            pt_orig = pt_orig(1:2) / pt_orig(3);

            % Validate: point should be within reasonable image bounds
            if pt_orig(1) > 0 && pt_orig(2) > 0 && pt_orig(1) < 5000 && pt_orig(2) < 5000
                nodal_points = [nodal_points; pt_orig', i, j];
                fprintf('  N(%d,%d) at (%.1f, %.1f)\n', i, j, pt_orig(1), pt_orig(2));
            end
        end
    end
end

fprintf('Found %d nodal points\n', size(nodal_points, 1));
end

%% ========== HELPER FUNCTION ==========

function [pt, found] = find_polyline_intersection(pts1, pts2)
% FIND_POLYLINE_INTERSECTION - Find where two polylines cross
%   Uses strict segment intersection (t1, t2 both in [0,1])
%
%   This is the geometrically correct method:
%   - No curve fitting required
%   - Finds actual crossing point of discrete arc samples
%   - Works directly with the user-selected arc points

pt = [];
found = false;

n1 = size(pts1, 1);
n2 = size(pts2, 1);

if n1 < 2 || n2 < 2
    return;
end

% Check each segment of arc1 against each segment of arc2
for i = 1:n1-1
    p1a = pts1(i, :);
    p1b = pts1(i+1, :);

    for j = 1:n2-1
        p2a = pts2(j, :);
        p2b = pts2(j+1, :);

        % Line segment intersection
        d1 = p1b - p1a;
        d2 = p2b - p2a;
        d12 = p1a - p2a;

        denom = d1(1)*d2(2) - d1(2)*d2(1);

        if abs(denom) < 1e-10
            continue;  % Parallel segments
        end

        t1 = (d2(1)*d12(2) - d2(2)*d12(1)) / denom;
        t2 = (d1(1)*d12(2) - d1(2)*d12(1)) / denom;

        % STRICT: intersection must be within BOTH segments
        if t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1
            pt = p1a + t1 * d1;
            found = true;
            return;  % Found true crossing
        end
    end
end
end

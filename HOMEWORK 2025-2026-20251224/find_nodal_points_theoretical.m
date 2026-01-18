function nodal_points = find_nodal_points_theoretical(arcs_A, arcs_B, lines_trans, H_R, v_axis)
% FIND_NODAL_POINTS_THEORETICAL - Find nodal points using symmetry plane theory
%
% Theory:
%   - Symmetry planes perpendicular to cylinder axis become HORIZONTAL in rectified space
%   - Each transversal line defines the Y-coordinate of a symmetry plane in rectified space
%   - Nodal point N_ij: where arc A_i and B_j cross the SAME horizontal level in rectified
%
% Inputs:
%   arcs_A, arcs_B: Cell arrays of arc points (image coords)
%   lines_trans: Transversal lines - define Y-levels of symmetry planes
%   H_R: Rectification homography
%   v_axis: Vanishing point of cylinder axis

nodal_points = [];

if isempty(arcs_A) || isempty(arcs_B)
    fprintf('No arcs provided.\n');
    return;
end

if nargin < 4 || isempty(H_R)
    H_R = eye(3);
end
H_R_inv = inv(H_R);

fprintf('Finding nodal points using symmetry planes in rectified space...\n');

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

%% Get Y-levels from transversal lines (symmetry plane heights in rectified space)
% In rectified space, transversal lines should be HORIZONTAL
% Their Y-coordinate defines the symmetry plane

y_levels = [];
if ~isempty(lines_trans)
    for t = 1:length(lines_trans)
        L = lines_trans(t);
        % Transform transversal midpoint to rectified space
        mid = (L.p1 + L.p2) / 2;
        mid_rect = H_R * [mid, 1]';
        mid_rect = mid_rect(1:2) / mid_rect(3);
        y_levels = [y_levels; mid_rect(2), t];
    end
    fprintf('  %d symmetry planes at Y = %.1f to %.1f (rectified)\n', ...
        length(lines_trans), min(y_levels(:,1)), max(y_levels(:,1)));
end

%% Find arc crossings at each Y-level (symmetry plane)
Y_TOLERANCE = 50;  % pixels in rectified space - how close Y must be

for lev = 1:size(y_levels, 1)
    y_plane = y_levels(lev, 1);
    trans_idx = y_levels(lev, 2);

    % Find where each A-arc crosses this Y-level
    crossings_A = [];  % [arc_idx, x, y]
    for i = 1:length(arcs_A_rect)
        pts = arcs_A_rect{i};
        [x_cross, found] = find_y_crossing(pts, y_plane);
        if found
            crossings_A = [crossings_A; i, x_cross, y_plane];
        end
    end

    % Find where each B-arc crosses this Y-level
    crossings_B = [];  % [arc_idx, x, y]
    for j = 1:length(arcs_B_rect)
        pts = arcs_B_rect{j};
        [x_cross, found] = find_y_crossing(pts, y_plane);
        if found
            crossings_B = [crossings_B; j, x_cross, y_plane];
        end
    end

    fprintf('  Y=%.0f (Trans%d): %d A, %d B crossings\n', y_plane, trans_idx, ...
        size(crossings_A,1), size(crossings_B,1));

    % Match crossings - where A and B cross at similar X position
    X_PROXIMITY = 20;  % pixels

    for a = 1:size(crossings_A, 1)
        x_A = crossings_A(a, 2);
        arc_A_idx = crossings_A(a, 1);

        for b = 1:size(crossings_B, 1)
            x_B = crossings_B(b, 2);
            arc_B_idx = crossings_B(b, 1);

            x_dist = abs(x_A - x_B);

            if x_dist < X_PROXIMITY
                % Found nodal point in rectified space
                nodal_rect = [(x_A + x_B)/2, y_plane];

                % Transform back to original image
                nodal_orig = H_R_inv * [nodal_rect, 1]';
                nodal_orig = nodal_orig(1:2) / nodal_orig(3);

                if nodal_orig(1) > 0 && nodal_orig(2) > 0
                    nodal_points = [nodal_points; nodal_orig', arc_A_idx, arc_B_idx];
                    fprintf('    N(%d,%d) at (%.1f, %.1f) [x_dist=%.1f]\n', ...
                        arc_A_idx, arc_B_idx, nodal_orig(1), nodal_orig(2), x_dist);
                end
            end
        end
    end
end

%% Also check direct arc intersections (may catch nodes between defined Y-levels)
fprintf('  Checking direct arc intersections...\n');
for i = 1:length(arcs_A_rect)
    for j = 1:length(arcs_B_rect)
        [pt_rect, found] = polyline_intersection(arcs_A_rect{i}, arcs_B_rect{j});
        if found
            pt_orig = H_R_inv * [pt_rect(:); 1];
            pt_orig = pt_orig(1:2) / pt_orig(3);
            if pt_orig(1) > 0 && pt_orig(2) > 0
                % Check if not already found
                is_duplicate = false;
                for k = 1:size(nodal_points, 1)
                    if nodal_points(k,3) == i && nodal_points(k,4) == j
                        is_duplicate = true;
                        break;
                    end
                end
                if ~is_duplicate
                    nodal_points = [nodal_points; pt_orig', i, j];
                    fprintf('    N(%d,%d) at (%.1f, %.1f) [direct]\n', i, j, pt_orig(1), pt_orig(2));
                end
            end
        end
    end
end

fprintf('Found %d nodal points\n', size(nodal_points, 1));
end

%% ========== HELPER FUNCTIONS ==========

function [x_cross, found] = find_y_crossing(pts, y_target)
% FIND_Y_CROSSING - Find where a polyline crosses a horizontal line Y=y_target

x_cross = [];
found = false;

for k = 1:size(pts, 1)-1
    y1 = pts(k, 2);
    y2 = pts(k+1, 2);

    % Check if segment crosses the Y-level
    if (y1 <= y_target && y2 >= y_target) || (y1 >= y_target && y2 <= y_target)
        if abs(y2 - y1) > 1e-6
            t = (y_target - y1) / (y2 - y1);
            if t >= 0 && t <= 1
                x1 = pts(k, 1);
                x2 = pts(k+1, 1);
                x_cross = x1 + t * (x2 - x1);
                found = true;
                return;
            end
        end
    end
end
end

function [pt, found] = polyline_intersection(pts1, pts2)
pt = [];
found = false;

for i = 1:size(pts1, 1)-1
    p1a = pts1(i, :);
    p1b = pts1(i+1, :);
    for j = 1:size(pts2, 1)-1
        p2a = pts2(j, :);
        p2b = pts2(j+1, :);

        d1 = p1b - p1a;
        d2 = p2b - p2a;
        denom = d1(1)*d2(2) - d1(2)*d2(1);
        if abs(denom) < 1e-10, continue; end

        d12 = p1a - p2a;
        t1 = (d2(1)*d12(2) - d2(2)*d12(1)) / denom;
        t2 = (d1(1)*d12(2) - d1(2)*d12(1)) / denom;

        if t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1
            pt = p1a + t1 * d1;
            found = true;
            return;
        end
    end
end
end

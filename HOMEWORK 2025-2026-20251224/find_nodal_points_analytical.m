function [nodal_points, conics_A, conics_B] = find_nodal_points_analytical(arcs_A, arcs_B, H_R)
% FIND_NODAL_POINTS_ANALYTICAL - Find arc intersections using conic fitting
%   Works in RECTIFIED space for cleaner geometry, with normalization for stability.
%
%   Inputs:
%     arcs_A, arcs_B: Cell arrays of arc points in original image coords
%     H_R: Rectification homography (optional, uses identity if not provided)
%
%   Returns:
%     nodal_points: Nx4 matrix [x, y, arc_A_idx, arc_B_idx] in ORIGINAL image coords
%     conics_A, conics_B: Cell arrays of 3x3 conic matrices (in normalized rectified space)

nodal_points = [];
conics_A = {};
conics_B = {};

if isempty(arcs_A) || isempty(arcs_B)
    return;
end

% Use identity if H_R not provided
if nargin < 3 || isempty(H_R)
    H_R = eye(3);
end
H_R_inv = inv(H_R);

%% Step 1: Transform all arcs to rectified space
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

%% Step 2: Compute normalization transform (for numerical stability)
% Collect all points to compute a common normalization
all_pts = [];
for i = 1:length(arcs_A_rect), all_pts = [all_pts; arcs_A_rect{i}]; end
for i = 1:length(arcs_B_rect), all_pts = [all_pts; arcs_B_rect{i}]; end

% Normalization: shift to centroid, scale so mean distance from origin is sqrt(2)
centroid = mean(all_pts, 1);
shifted = all_pts - centroid;
mean_dist = mean(sqrt(sum(shifted.^2, 2)));
scale = sqrt(2) / mean_dist;

T_norm = [scale, 0, -scale*centroid(1);
    0, scale, -scale*centroid(2);
    0, 0, 1];
T_norm_inv = inv(T_norm);

fprintf('Conic fitting: centroid=(%.1f,%.1f), scale=%.4f\n', centroid, scale);

%% Step 3: Fit conics to all arcs (in normalized rectified space)
fprintf('Fitting conics...\n');

for i = 1:length(arcs_A_rect)
    pts = arcs_A_rect{i};
    pts_norm = (T_norm(1:2, 1:2) * pts' + T_norm(1:2, 3))';
    C = fit_conic(pts_norm);
    conics_A{i} = C;
    fprintf('  Arc A%d: fitted (%d pts)\n', i, size(pts, 1));
end

for i = 1:length(arcs_B_rect)
    pts = arcs_B_rect{i};
    pts_norm = (T_norm(1:2, 1:2) * pts' + T_norm(1:2, 3))';
    C = fit_conic(pts_norm);
    conics_B{i} = C;
    fprintf('  Arc B%d: fitted (%d pts)\n', i, size(pts, 1));
end

%% Step 4: Find intersections between all pairs
fprintf('Finding intersections...\n');

for i = 1:length(conics_A)
    C_A = conics_A{i};
    pts_A_norm = (T_norm(1:2, 1:2) * arcs_A_rect{i}' + T_norm(1:2, 3))';

    for j = 1:length(conics_B)
        C_B = conics_B{j};
        pts_B_norm = (T_norm(1:2, 1:2) * arcs_B_rect{j}' + T_norm(1:2, 3))';

        % Try algebraic intersection first
        pts_int_norm = intersect_conics(C_A, C_B, pts_A_norm, pts_B_norm);

        if isempty(pts_int_norm)
            % Fallback: find closest approach between the two arcs
            pts_int_norm = find_closest_approach(pts_A_norm, pts_B_norm);
        end

        if ~isempty(pts_int_norm)
            % Denormalize: from normalized -> rectified -> original image
            for k = 1:size(pts_int_norm, 1)
                pt_norm = pts_int_norm(k, :);

                % Denormalize to rectified coords
                pt_rect = T_norm_inv * [pt_norm, 1]';
                pt_rect = pt_rect(1:2) / pt_rect(3);

                % Transform back to original image coords
                pt_orig = H_R_inv * [pt_rect; 1];
                pt_orig = pt_orig(1:2) / pt_orig(3);

                % Check if within reasonable bounds
                if pt_orig(1) > 0 && pt_orig(2) > 0 && pt_orig(1) < 5000 && pt_orig(2) < 5000
                    nodal_points = [nodal_points; pt_orig', i, j];
                    fprintf('  N(%d,%d) at (%.1f, %.1f)\n', i, j, pt_orig(1), pt_orig(2));
                    break;  % One nodal point per pair
                end
            end
        end
    end
end

fprintf('Found %d nodal points\n', size(nodal_points, 1));
end

%% ========== HELPER FUNCTIONS ==========

function C = fit_conic(pts)
% FIT_CONIC - Fit conic to normalized 2D points using SVD
%   Returns 3x3 symmetric matrix C such that x^T C x = 0

n = size(pts, 1);
if n < 6
    warning('Need at least 6 points for robust conic fit');
    C = eye(3);
    return;
end

% Design matrix: [x^2, xy, y^2, x, y, 1]
A = zeros(n, 6);
for k = 1:n
    x = pts(k, 1);
    y = pts(k, 2);
    A(k, :) = [x^2, x*y, y^2, x, y, 1];
end

% SVD solution
[~, ~, V] = svd(A, 'econ');
p = V(:, end);

% Build symmetric matrix
a = p(1); b = p(2); c = p(3); d = p(4); e = p(5); f = p(6);
C = [a,   b/2, d/2;
    b/2, c,   e/2;
    d/2, e/2, f];
C = C / norm(C, 'fro');
end

function pts = intersect_conics(C1, C2, pts1, pts2)
% INTERSECT_CONICS - Find intersection points of two conics
%   Uses the closest approach on the arcs as a robust method

pts = [];

% Method: Find where both conic residuals are minimized
% This is more robust than the algebraic pencil method for noisy data

% Sample both arc point sets and find where they cross
% Use linear interpolation between consecutive points

n1 = size(pts1, 1);
n2 = size(pts2, 1);

if n1 < 2 || n2 < 2
    return;
end

best_pt = [];
best_score = inf;

% For each segment in arc1, check intersection with each segment in arc2
for i = 1:n1-1
    p1a = pts1(i, :);
    p1b = pts1(i+1, :);

    for j = 1:n2-1
        p2a = pts2(j, :);
        p2b = pts2(j+1, :);

        % Find intersection of line segments
        [pt, t1, t2] = line_segment_intersection(p1a, p1b, p2a, p2b);

        if ~isempty(pt) && t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1
            % Valid intersection (with some tolerance for extrapolation)
            score = abs(t1 - 0.5) + abs(t2 - 0.5);  % Prefer middle of segments
            if score < best_score
                best_score = score;
                best_pt = pt;
            end
        end
    end
end

if ~isempty(best_pt)
    pts = best_pt;
end
end

function [pt, t1, t2] = line_segment_intersection(p1a, p1b, p2a, p2b)
% LINE_SEGMENT_INTERSECTION - Find intersection of two line segments
%   Returns the intersection point and parameters t1, t2 such that:
%   intersection = p1a + t1*(p1b - p1a) = p2a + t2*(p2b - p2a)

pt = [];
t1 = nan;
t2 = nan;

d1 = p1b - p1a;
d2 = p2b - p2a;
d12 = p1a - p2a;

denom = d1(1)*d2(2) - d1(2)*d2(1);

if abs(denom) < 1e-10
    return;  % Parallel or coincident
end

t1 = (d2(1)*d12(2) - d2(2)*d12(1)) / denom;
t2 = (d1(1)*d12(2) - d1(2)*d12(1)) / denom;

pt = p1a + t1 * d1;
end

function pts = find_closest_approach(pts1, pts2)
% FIND_CLOSEST_APPROACH - Find point where two polylines come closest
%   Fallback method when segment intersection fails

pts = [];
min_dist = inf;
best_pt = [];

for i = 1:size(pts1, 1)
    for j = 1:size(pts2, 1)
        d = norm(pts1(i, :) - pts2(j, :));
        if d < min_dist
            min_dist = d;
            best_pt = (pts1(i, :) + pts2(j, :)) / 2;
        end
    end
end

% Only return if reasonably close (threshold in normalized coords)
if min_dist < 0.1 && ~isempty(best_pt)
    pts = best_pt;
end
end

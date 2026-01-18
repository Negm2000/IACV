% Submission.m - IACV Homework 2025-2026
% 3D Reconstruction of a Cylindric Vault
% Student submission - single file solution

clear;
close all;
clc;

%% ========== CONFIGURATION ==========
imageFile = 'San Maurizio.jpg';
featureFile = 'features_fixed.mat';

%% ========== 1. LOAD IMAGE AND FEATURES ==========
fprintf('=== Loading ===\n');
img = imread(imageFile);
[rows, cols, ~] = size(img);
load(featureFile);
fprintf('Loaded: %d vertical lines, %d axis lines, %d transversal lines\n', ...
    length(lines_v), length(lines_axis), length(lines_trans));
fprintf('Loaded: %d arcs A, %d arcs B\n', length(arcs_A), length(arcs_B));

%% ========== 2. COMPUTE VANISHING POINTS ==========
fprintf('\n=== Computing Vanishing Points ===\n');

% Image center for normalization
cx = cols / 2;
cy = rows / 2;
scale = max(cx, cy);

% Normalization matrix
T_norm_to_pixel = [scale, 0, cx; 0, scale, cy; 0, 0, 1];

% Function to convert two points to homogeneous line (normalized coordinates)
% Line equation: l = p1 x p2 (cross product)
get_normalized_line = @(p1, p2) cross(...
    [(p1(1)-cx)/scale, (p1(2)-cy)/scale, 1], ...
    [(p2(1)-cx)/scale, (p2(2)-cy)/scale, 1])';

% --- Vertical vanishing point ---
if length(lines_v) >= 2
    A_vert = zeros(length(lines_v), 3);
    for i = 1:length(lines_v)
        line_hom = get_normalized_line(lines_v(i).p1, lines_v(i).p2);
        A_vert(i, :) = line_hom';
    end
    [~, ~, V] = svd(A_vert);
    v_vert_normalized = V(:, end)';
    v_vert_normalized = v_vert_normalized / v_vert_normalized(3);
else
    error('Need at least 2 vertical lines');
end

% --- Axis vanishing point (along the barrel vault) ---
% Combine axis lines and apical line
all_axis_lines = {};
for i = 1:length(lines_axis)
    all_axis_lines{end+1} = get_normalized_line(lines_axis(i).p1, lines_axis(i).p2);
end
if ~isempty(line_apical)
    all_axis_lines{end+1} = get_normalized_line(line_apical.p1, line_apical.p2);
end

A_axis = zeros(length(all_axis_lines), 3);
for i = 1:length(all_axis_lines)
    A_axis(i, :) = all_axis_lines{i}';
end
[~, ~, V] = svd(A_axis);
v_axis_normalized = V(:, end)';
v_axis_normalized = v_axis_normalized / v_axis_normalized(3);

% --- Transversal vanishing point ---
if length(lines_trans) >= 2
    A_trans = zeros(length(lines_trans), 3);
    for i = 1:length(lines_trans)
        line_hom = get_normalized_line(lines_trans(i).p1, lines_trans(i).p2);
        A_trans(i, :) = line_hom';
    end
    [~, ~, V] = svd(A_trans);
    v_trans_normalized = V(:, end)';
    v_trans_normalized = v_trans_normalized / v_trans_normalized(3);
else
    error('Need at least 2 transversal lines');
end

% Convert to pixel coordinates
v_vert_pixel = (T_norm_to_pixel * v_vert_normalized(:))';
v_vert_pixel = v_vert_pixel / v_vert_pixel(3);

v_axis_pixel = (T_norm_to_pixel * v_axis_normalized(:))';
v_axis_pixel = v_axis_pixel / v_axis_pixel(3);

v_trans_pixel = (T_norm_to_pixel * v_trans_normalized(:))';
v_trans_pixel = v_trans_pixel / v_trans_pixel(3);

% Vanishing line (line at infinity in the plane perpendicular to vertical)
vanishing_line = cross(v_vert_pixel, v_trans_pixel);
vanishing_line = vanishing_line / norm(vanishing_line(1:2));

fprintf('Vertical VP: (%.1f, %.1f)\n', v_vert_pixel(1)/v_vert_pixel(3), v_vert_pixel(2)/v_vert_pixel(3));
fprintf('Axis VP: (%.1f, %.1f)\n', v_axis_pixel(1)/v_axis_pixel(3), v_axis_pixel(2)/v_axis_pixel(3));
fprintf('Trans VP: (%.1f, %.1f)\n', v_trans_pixel(1)/v_trans_pixel(3), v_trans_pixel(2)/v_trans_pixel(3));

%% ========== 3. METRIC RECTIFICATION ==========
fprintf('\n=== Metric Rectification ===\n');

% Step 1: Affine rectification - map vanishing line to infinity
l = vanishing_line(:) / vanishing_line(3);
H_affine = [1, 0, 0; 0, 1, 0; l(1), l(2), 1];

% Step 2: Metric rectification using perpendicular directions
v_vert_affine = H_affine * v_vert_pixel(:);
v_trans_affine = H_affine * v_trans_pixel(:);

dir_vert = v_vert_affine(1:2) / norm(v_vert_affine(1:2));
dir_trans = v_trans_affine(1:2) / norm(v_trans_affine(1:2));

M = [dir_vert, dir_trans];
S_metric = [0, 1; 1, 0] * inv(M);
H_metric = [S_metric, [0; 0]; 0, 0, 1];

H_rectify = H_metric * H_affine;

% Compute output bounds by sampling grid points
[xx, yy] = meshgrid(linspace(1, cols, 40), linspace(1, rows, 40));
grid_points = [xx(:), yy(:), ones(numel(xx), 1)]';

grid_vals = H_affine(3, :) * grid_points;
side = sign(mean(grid_vals));
if side == 0
    side = 1;
end

% Filter points on correct side of vanishing line
margin = 0.15 * abs(max(grid_vals) - min(grid_vals));
valid_mask = (sign(grid_vals) == side) & (abs(grid_vals) > margin);
if sum(valid_mask) < 10
    valid_mask = (sign(grid_vals) == side);
end

pts_transformed = H_rectify * grid_points(:, valid_mask);
pts_transformed = pts_transformed(1:2, :) ./ pts_transformed(3, :);

% Bounding box
min_x = min(pts_transformed(1, :));
max_x = max(pts_transformed(1, :));
min_y = min(pts_transformed(2, :));
max_y = max(pts_transformed(2, :));

rect_width = max_x - min_x;
rect_height = max_y - min_y;

% Limit aspect ratio
if rect_height > 5 * rect_width
    rect_height = 5 * rect_width;
end
if rect_width > 5 * rect_height
    rect_width = 5 * rect_height;
end

% Scale to fit target width
target_width = 2000;
scale_factor = target_width / rect_width;

% Translation and scale matrix
T_output = [scale_factor, 0, -scale_factor*min_x + 1; ...
    0, scale_factor, -scale_factor*min_y + 1; ...
    0, 0, 1];

H_final = T_output * H_rectify;

output_size = [ceil(scale_factor * rect_height), target_width];
img_rectified = imwarp(img, projective2d(H_final'), 'OutputView', imref2d(output_size));

fprintf('Rectified image size: %d x %d\n', output_size(2), output_size(1));

% Display a sample non-apical nodal point (will be computed after finding intersections)
% This is shown later after nodal points are found

%% ========== 4. CAMERA CALIBRATION ==========
fprintf('\n=== Camera Calibration ===\n');

% Using 3 orthogonal vanishing points to compute K
% Theory: vi^T * omega * vj = 0 for perpendicular directions
% omega = K^(-T) * K^(-1) is the Image of Absolute Conic

v1 = v_vert_normalized(:);
v2 = v_axis_normalized(:);
v3 = v_trans_normalized(:);

% Build system of equations
% omega parametrization (zero skew): [w1, 0, w4; 0, w3, w5; w4, w5, 1]
% Constraint: v1^T * omega * v2 = 0 gives linear equation in w1, w3, w4, w5

A_calib = zeros(4, 4);
b_calib = zeros(4, 1);

% v1 perpendicular to v2
A_calib(1, :) = [v1(1)*v2(1), v1(2)*v2(2), v1(1)*v2(3)+v1(3)*v2(1), v1(2)*v2(3)+v1(3)*v2(2)];
b_calib(1) = -v1(3)*v2(3);

% v1 perpendicular to v3
A_calib(2, :) = [v1(1)*v3(1), v1(2)*v3(2), v1(1)*v3(3)+v1(3)*v3(1), v1(2)*v3(3)+v1(3)*v3(2)];
b_calib(2) = -v1(3)*v3(3);

% v2 perpendicular to v3
A_calib(3, :) = [v2(1)*v3(1), v2(2)*v3(2), v2(1)*v3(3)+v2(3)*v3(1), v2(2)*v3(3)+v2(3)*v3(2)];
b_calib(3) = -v2(3)*v3(3);

% Square pixels assumption: w1 = w3
A_calib(4, :) = [1, -1, 0, 0];
b_calib(4) = 0;

% Solve for omega parameters
x = A_calib \ b_calib;
w1 = x(1);
w3 = x(2);
w4 = x(3);
w5 = x(4);

omega = [w1, 0, w4; 0, w3, w5; w4, w5, 1];

% Extract K via Cholesky decomposition
try
    L = chol(omega, 'lower');
    K_normalized = inv(L');
    K_normalized = K_normalized / K_normalized(3, 3);
catch
    % Handle non-positive-definite case
    [U, D, ~] = svd(omega);
    D_positive = abs(D);
    omega_fixed = U * D_positive * U';
    L = chol(omega_fixed, 'lower');
    K_normalized = inv(L');
    K_normalized = K_normalized / K_normalized(3, 3);
end

% Denormalize to pixel coordinates
K = T_norm_to_pixel * K_normalized;
K = K / K(3, 3);

% Ensure positive focal lengths
if K(1, 1) < 0
    K = -K;
    K = K / K(3, 3);
end

fprintf('Focal length: fx = %.1f, fy = %.1f\n', K(1,1), K(2,2));
fprintf('Principal point: (%.1f, %.1f)\n', K(1,3), K(2,3));

%% ========== 5. FIND NODAL POINTS ==========
fprintf('\n=== Finding Nodal Points ===\n');

H_inv = inv(H_final);
nodal_points = [];

% Find intersections between all pairs of arcs
for i = 1:length(arcs_A)
    % Transform arc A to rectified space
    pts_A = arcs_A{i};
    pts_A_hom = [pts_A, ones(size(pts_A, 1), 1)]';
    pts_A_rect = H_final * pts_A_hom;
    pts_A_rect = (pts_A_rect(1:2, :) ./ pts_A_rect(3, :))';

    for j = 1:length(arcs_B)
        % Transform arc B to rectified space
        pts_B = arcs_B{j};
        pts_B_hom = [pts_B, ones(size(pts_B, 1), 1)]';
        pts_B_rect = H_final * pts_B_hom;
        pts_B_rect = (pts_B_rect(1:2, :) ./ pts_B_rect(3, :))';

        % Find intersection
        intersection_pt = find_intersection(pts_A_rect, pts_B_rect);

        if ~isempty(intersection_pt)
            % Transform back to original image
            pt_orig = H_inv * [intersection_pt(:); 1];
            pt_orig = pt_orig(1:2) / pt_orig(3);

            if pt_orig(1) > 0 && pt_orig(2) > 0 && pt_orig(1) < 10000
                nodal_points = [nodal_points; pt_orig', i, j];
            end
        end
    end
end

fprintf('Found %d intersections\n', size(nodal_points, 1));

% Detect index shift using apical line
shift = 0;
if ~isempty(line_apical) && ~isempty(nodal_points)
    p1 = line_apical.p1;
    p2 = line_apical.p2;

    % Line equation: ax + by + c = 0
    a = p1(2) - p2(2);
    b = p2(1) - p1(1);
    c = -a*p1(1) - b*p1(2);
    line_norm = sqrt(a^2 + b^2);

    % Try different shifts and find best alignment with apical line
    best_count = 0;
    best_shift = 0;

    for k = -5:5
        % Get points where B_index - A_index = k
        candidates = nodal_points(nodal_points(:,4) - nodal_points(:,3) == k, :);

        if isempty(candidates)
            continue;
        end

        % Count points close to apical line
        count = 0;
        for idx = 1:size(candidates, 1)
            pt = candidates(idx, 1:2);
            distance = abs(a*pt(1) + b*pt(2) + c) / line_norm;
            if distance < 50  % 50 pixel tolerance
                count = count + 1;
            end
        end

        if count > best_count
            best_count = count;
            best_shift = k;
        end
    end

    shift = best_shift;
    fprintf('Detected index shift: %d\n', shift);

    % Apply shift correction
    nodal_points(:, 4) = nodal_points(:, 4) - shift;
end

% Print apical nodes (where A_index == B_index)
fprintf('\n--- Apical Nodal Points (on cylinder apex) ---\n');
for i = 1:size(nodal_points, 1)
    if nodal_points(i, 3) == nodal_points(i, 4)
        fprintf('APICAL Node N(%d,%d) at (%.0f, %.0f)\n', ...
            nodal_points(i, 3), nodal_points(i, 4), nodal_points(i, 1), nodal_points(i, 2));
    end
end

% Print a non-apical nodal point (requirement 3: compute position of non-apical nodal point)
% Pick an interior nodal point (both indices > 0) to avoid edge outliers
fprintf('\n--- Non-Apical Nodal Point (on cylinder surface) ---\n');
non_apical_mask = (nodal_points(:, 3) ~= nodal_points(:, 4)) & ...
    (nodal_points(:, 3) > 0) & (nodal_points(:, 4) > 0);
non_apical_idx = find(non_apical_mask, 1);
if ~isempty(non_apical_idx)
    np = nodal_points(non_apical_idx, :);
    fprintf('Non-apical Node N(%d,%d) at image coords: (%.1f, %.1f)\n', np(3), np(4), np(1), np(2));
    % Transform to rectified space
    np_rect = H_final * [np(1:2), 1]';
    np_rect = np_rect(1:2) / np_rect(3);
    fprintf('Same point in rectified plane: (%.1f, %.1f)\n', np_rect(1), np_rect(2));
end

%% ========== 6. 3D RECONSTRUCTION ==========
fprintf('\n=== 3D Reconstruction ===\n');

K_inv = inv(K);

% Compute 3D directions from vanishing points
vert_direction = K_inv * v_vert_pixel(:);
vert_direction = vert_direction / norm(vert_direction);

axis_direction = K_inv * v_axis_pixel(:);
axis_direction = axis_direction / norm(axis_direction);
axis_direction = axis_direction(:)';  % Make row vector

trans_direction = K_inv * v_trans_pixel(:);
trans_direction = trans_direction / norm(trans_direction);

fprintf('3D Direction vectors (in camera frame):\n');
fprintf('  Vertical:    [%.4f, %.4f, %.4f]\n', vert_direction);
fprintf('  Axis:        [%.4f, %.4f, %.4f]\n', axis_direction);
fprintf('  Transversal: [%.4f, %.4f, %.4f]\n', trans_direction);

% Arc spacing constraint (known distance between arcs)
d = 1;

% Get apical nodes (where A_index == B_index)
apical_mask = nodal_points(:, 3) == nodal_points(:, 4);
apical_nodes = nodal_points(apical_mask, :);
non_apical_nodes = nodal_points(~apical_mask, :);

% Sort apical nodes by index
[~, sort_idx] = sort(apical_nodes(:, 3));
apical_nodes = apical_nodes(sort_idx, :);

% Compute reference depth using two apical nodes
n1_2d = apical_nodes(1, 1:2);
n2_2d = apical_nodes(2, 1:2);

ray1 = K_inv * [n1_2d, 1]';
ray1 = ray1 / norm(ray1);

ray2 = K_inv * [n2_2d, 1]';
ray2 = ray2 / norm(ray2);

a1 = dot(ray1, axis_direction);
a2 = dot(ray2, axis_direction);

idx1 = apical_nodes(1, 3);
idx2 = apical_nodes(2, 3);
delta_d = (idx2 - idx1) * d;

% Solve for reference depth
lambda_ref = abs(delta_d / (a2 - a1));
fprintf('Reference depth (lambda): %.3f\n', lambda_ref);

% 3D position of apex
P_apex = lambda_ref * ray1(:)';

% Compute cylinder radius using all non-apical nodes
num_nodes = size(non_apical_nodes, 1);
lambdas = zeros(num_nodes, 1);

for i = 1:num_nodes
    node_2d = non_apical_nodes(i, 1:2);
    ray_i = K_inv * [node_2d, 1]';
    ray_i = ray_i / norm(ray_i);
    a_i = dot(ray_i, axis_direction);

    if abs(a_i) > 1e-6
        lambdas(i) = lambda_ref * a1 / a_i;
    else
        lambdas(i) = Inf;
    end
end

% Outlier detection using median
finite_lambdas = lambdas(isfinite(lambdas));
median_lambda = median(finite_lambdas);
mad_lambda = median(abs(finite_lambdas - median_lambda));
threshold = max(3 * mad_lambda, 0.5 * lambda_ref);

inlier_mask = abs(lambdas - median_lambda) < threshold & isfinite(lambdas);

% Compute radius from all inlier nodes
radius_estimates = [];
for i = 1:num_nodes
    if inlier_mask(i)
        node_2d = non_apical_nodes(i, 1:2);
        ray_i = K_inv * [node_2d, 1]';
        ray_i = ray_i / norm(ray_i);

        P_side = lambdas(i) * ray_i(:)';

        diff = P_side - P_apex;
        diff_perp = diff - dot(diff, axis_direction) * axis_direction;

        radius_estimates = [radius_estimates; norm(diff_perp)];
    end
end

R_cylinder = median(radius_estimates);
fprintf('Cylinder radius R: %.3f\n', R_cylinder);

% Cylinder axis position (below apex)
P_axis = P_apex - R_cylinder * vert_direction(:)';

% Print cylinder axis localization (requirement 6)
fprintf('\n--- Cylinder Axis Localization (wrt Camera) ---\n');
fprintf('Axis passes through point: [%.4f, %.4f, %.4f]\n', P_axis);
fprintf('Axis direction vector:     [%.4f, %.4f, %.4f]\n', axis_direction);

% Reconstruct all arc points
points_3D = {};

for i = 1:length(arcs_A)
    arc_2d = arcs_A{i};
    pts_3d = [];

    for j = 1:size(arc_2d, 1)
        ray = K_inv * [arc_2d(j, :), 1]';
        ray = ray(:)' / norm(ray);

        lambda = intersect_ray_cylinder(ray, P_axis, axis_direction, R_cylinder, lambda_ref);

        if lambda > 0 && isfinite(lambda)
            pt_3d = lambda * ray;
            pts_3d = [pts_3d; pt_3d];
        end
    end

    if ~isempty(pts_3d)
        points_3D{end+1} = struct('pts', pts_3d, 'arc', 'A');
    end
end

for i = 1:length(arcs_B)
    arc_2d = arcs_B{i};
    pts_3d = [];

    for j = 1:size(arc_2d, 1)
        ray = K_inv * [arc_2d(j, :), 1]';
        ray = ray(:)' / norm(ray);

        lambda = intersect_ray_cylinder(ray, P_axis, axis_direction, R_cylinder, lambda_ref);

        if lambda > 0 && isfinite(lambda)
            pt_3d = lambda * ray;
            pts_3d = [pts_3d; pt_3d];
        end
    end

    if ~isempty(pts_3d)
        points_3D{end+1} = struct('pts', pts_3d, 'arc', 'B');
    end
end

% Verify reconstruction
centroid1 = mean(points_3D{1}.pts, 1);
centroid2 = mean(points_3D{2}.pts, 1);
measured_spacing = abs(dot(centroid2 - centroid1, axis_direction));

fprintf('\n%d arcs reconstructed\n', length(points_3D));
fprintf('Arc spacing verification: %.3f (target: 1.0)\n', measured_spacing);

% Print 3D coordinates of a dozen points from one arc (requirement 5)
fprintf('\n--- 3D Coordinates of 12 Points from Arc A1 (diagonal arc) ---\n');
arc1_pts = points_3D{1}.pts;
num_to_show = min(12, size(arc1_pts, 1));
fprintf('Point |     X     |     Y     |     Z     \n');
fprintf('------+-----------+-----------+-----------\n');
for i = 1:num_to_show
    fprintf('  %2d  | %9.4f | %9.4f | %9.4f\n', i, arc1_pts(i, 1), arc1_pts(i, 2), arc1_pts(i, 3));
end

%% ========== 7. VISUALIZATION ==========

% Figure 1: Original image with features
figure(1);
imshow(img);
title('Figure 1: Original Image with Extracted Features');
hold on;

% Draw vertical lines (black)
for i = 1:length(lines_v)
    L = lines_v(i);
    plot([L.p1(1), L.p2(1)], [L.p1(2), L.p2(2)], 'k-', 'LineWidth', 2);
end

% Draw axis lines (green)
for i = 1:length(lines_axis)
    L = lines_axis(i);
    plot([L.p1(1), L.p2(1)], [L.p1(2), L.p2(2)], 'g-', 'LineWidth', 2);
end

% Draw transversal lines (white)
for i = 1:length(lines_trans)
    L = lines_trans(i);
    plot([L.p1(1), L.p2(1)], [L.p1(2), L.p2(2)], 'w-', 'LineWidth', 2);
end

% Draw apical line (yellow)
if ~isempty(line_apical)
    plot([line_apical.p1(1), line_apical.p2(1)], ...
        [line_apical.p1(2), line_apical.p2(2)], 'y-', 'LineWidth', 3);
end

% Draw arcs
h_arcA = plot(arcs_A{1}(:,1), arcs_A{1}(:,2), 'c.-', 'MarkerSize', 8);
for i = 2:length(arcs_A)
    plot(arcs_A{i}(:,1), arcs_A{i}(:,2), 'c.-', 'MarkerSize', 8);
end
h_arcB = plot(arcs_B{1}(:,1), arcs_B{1}(:,2), 'm.-', 'MarkerSize', 8);
for i = 2:length(arcs_B)
    plot(arcs_B{i}(:,1), arcs_B{i}(:,2), 'm.-', 'MarkerSize', 8);
end

% Create proper legend with handles
h_vert = plot(NaN, NaN, 'k-', 'LineWidth', 2);
h_axis = plot(NaN, NaN, 'g-', 'LineWidth', 2);
h_trans = plot(NaN, NaN, 'w-', 'LineWidth', 2);
h_apic = plot(NaN, NaN, 'y-', 'LineWidth', 3);
lgd = legend([h_vert, h_axis, h_trans, h_apic, h_arcA, h_arcB], ...
    {'Vertical Lines', 'Axis Lines', 'Transversal Lines', 'Apical Line', 'Arcs A (cyan)', 'Arcs B (magenta)'}, ...
    'Location', 'best');
set(lgd, 'Color', [0.2 0.2 0.2], 'TextColor', 'w');  % Dark background with white text

% Figure 2: Vanishing Points and Vanishing Line
figure(2);
imshow(img);
title('Figure 2: Vanishing Points and Vanishing Line');
hold on;

% Draw the vanishing line
vl = vanishing_line;
if abs(vl(2)) > 1e-6
    x_range = [-cols, cols*2];
    y_range = (-vl(3) - vl(1)*x_range) / vl(2);
    plot(x_range, y_range, 'b-', 'LineWidth', 2);
end

% Collect all VP positions for axis bounds
all_vps = [];

% Draw lines converging to vertical VP
if abs(v_vert_pixel(3)) > 1e-6
    vp = v_vert_pixel(1:2) / v_vert_pixel(3);
    all_vps = [all_vps; vp];
    for i = 1:length(lines_v)
        plot([lines_v(i).p1(1), vp(1)], [lines_v(i).p1(2), vp(2)], 'k:', 'LineWidth', 0.5);
    end
    plot(vp(1), vp(2), 'ko', 'MarkerSize', 15, 'LineWidth', 3);
    text(vp(1)+20, vp(2), 'V_{vert}', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold');
end

% Draw lines converging to axis VP
if abs(v_axis_pixel(3)) > 1e-6
    vp = v_axis_pixel(1:2) / v_axis_pixel(3);
    all_vps = [all_vps; vp];
    for i = 1:length(lines_axis)
        plot([lines_axis(i).p1(1), vp(1)], [lines_axis(i).p1(2), vp(2)], 'g:', 'LineWidth', 0.5);
    end
    plot(vp(1), vp(2), 'go', 'MarkerSize', 15, 'LineWidth', 3);
    text(vp(1)+20, vp(2), 'V_{axis}', 'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold');
end

% Draw lines converging to transversal VP
if abs(v_trans_pixel(3)) > 1e-6
    vp = v_trans_pixel(1:2) / v_trans_pixel(3);
    all_vps = [all_vps; vp];
    for i = 1:length(lines_trans)
        plot([lines_trans(i).p1(1), vp(1)], [lines_trans(i).p1(2), vp(2)], 'r:', 'LineWidth', 0.5);
    end
    plot(vp(1), vp(2), 'ro', 'MarkerSize', 15, 'LineWidth', 3);
    text(vp(1)+20, vp(2), 'V_{trans}', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
end

% Extend axes to show vanishing points with margin
if ~isempty(all_vps)
    min_x = min([0; all_vps(:,1)]) - 200;
    max_x = max([cols; all_vps(:,1)]) + 200;
    min_y = min([0; all_vps(:,2)]) - 200;
    max_y = max([rows; all_vps(:,2)]) + 200;
    axis([min_x, max_x, min_y, max_y]);
end

% Add prominent grid and labels
axis on;
grid on;
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.6, 'GridColor', [0.2 0.2 0.2]);
set(gca, 'Layer', 'top'); % Put grid on top of the image
xlabel('Image X [pixels]');
ylabel('Image Y [pixels]');
box on;

% Legend
h_vl = plot(NaN, NaN, 'b-', 'LineWidth', 2);
h_vp_v = plot(NaN, NaN, 'ko', 'MarkerSize', 12, 'LineWidth', 2);
h_vp_a = plot(NaN, NaN, 'go', 'MarkerSize', 12, 'LineWidth', 2);
h_vp_t = plot(NaN, NaN, 'ro', 'MarkerSize', 12, 'LineWidth', 2);
legend([h_vl, h_vp_v, h_vp_a, h_vp_t], ...
    {'Vanishing Line', 'Vertical VP', 'Axis VP', 'Transversal VP'}, ...
    'Location', 'best');

% Figure 3: Rectified image with nodal points
figure(3);
imshow(img_rectified);
title('Figure 3: Metric Rectification with Nodal Points');
hold on;
% Show nodal points on rectified image
for i = 1:size(nodal_points, 1)
    pt_rect = H_final * [nodal_points(i, 1:2), 1]';
    pt_rect = pt_rect(1:2) / pt_rect(3);
    if nodal_points(i, 3) == nodal_points(i, 4)
        % Apical nodes - yellow
        plot(pt_rect(1), pt_rect(2), 'yo', 'MarkerSize', 12, 'LineWidth', 2);
        text(pt_rect(1)+10, pt_rect(2), sprintf('N%d%d', nodal_points(i,3), nodal_points(i,4)), ...
            'Color', 'y', 'FontSize', 10, 'FontWeight', 'bold');
    else
        % Non-apical nodes - red with labels
        plot(pt_rect(1), pt_rect(2), 'ro', 'MarkerSize', 8, 'LineWidth', 1);
        text(pt_rect(1)+10, pt_rect(2), sprintf('N%d%d', nodal_points(i,3), nodal_points(i,4)), ...
            'Color', 'r', 'FontSize', 8);
    end
end
% Add legend for nodal points
h_apical = plot(NaN, NaN, 'yo', 'MarkerSize', 12, 'LineWidth', 2);
h_nonapical = plot(NaN, NaN, 'ro', 'MarkerSize', 8, 'LineWidth', 1);
legend([h_apical, h_nonapical], {'Apical Nodes (N_ii)', 'Non-Apical Nodes (N_ij)'}, 'Location', 'best');

% Figure 4: Full 3D reconstruction
figure(4);
hold on;

% Plot all arcs
for i = 1:length(points_3D)
    pts = points_3D{i}.pts;
    if points_3D{i}.arc == 'A'
        color = 'c';
    else
        color = 'm';
    end
    plot3(pts(:,1), pts(:,2), pts(:,3), [color '.-'], 'LineWidth', 1.5, 'MarkerSize', 8);
end

% Plot cylinder axis
t = linspace(-3, 3, 30)';
axis_line = P_axis + t * axis_direction;
plot3(axis_line(:,1), axis_line(:,2), axis_line(:,3), 'g-', 'LineWidth', 3);

grid on;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title(sprintf('Figure 4: Full 3D Reconstruction (R = %.2f)', R_cylinder));
view(3);

% Figure 5: Multiple views of full reconstruction
figure(5);

views = {[0, 90], [0, 0], [90, 0], [45, 30]};
view_titles = {'Top View (XY)', 'Front View (XZ)', 'Side View (YZ)', 'Perspective'};

for v = 1:4
    subplot(2, 2, v);
    hold on;

    % Plot all arcs
    for i = 1:length(points_3D)
        pts = points_3D{i}.pts;
        if points_3D{i}.arc == 'A'
            color = 'c';
        else
            color = 'm';
        end
        plot3(pts(:,1), pts(:,2), pts(:,3), [color '.-']);
    end

    % Plot axis
    plot3(axis_line(:,1), axis_line(:,2), axis_line(:,3), 'g-', 'LineWidth', 2);

    view(views{v});
    axis equal;
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(view_titles{v});
end

sgtitle('Figure 5: 3D Reconstruction - Multiple Views');

% Figure 6: ONE curved arc with different views (Part 2 requirement 3)
figure(6);

single_arc = points_3D{1}.pts;  % First diagonal arc (Arc A1)

for v = 1:4
    subplot(2, 2, v);
    hold on;

    % Plot single arc
    plot3(single_arc(:,1), single_arc(:,2), single_arc(:,3), 'c.-', 'LineWidth', 2, 'MarkerSize', 10);

    % Plot cylinder axis for reference
    plot3(axis_line(:,1), axis_line(:,2), axis_line(:,3), 'g-', 'LineWidth', 2);

    view(views{v});
    axis equal;
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(view_titles{v});
end

sgtitle('Figure 6: Single Diagonal Arc (A1) - Different Views');

fprintf('\n=== Done ===\n');

%% ========== HELPER FUNCTIONS ==========

function pt = find_intersection(points1, points2)
% Find intersection between two polylines
pt = [];

n1 = size(points1, 1);
n2 = size(points2, 1);

if n1 < 2 || n2 < 2
    return;
end

for i = 1:n1-1
    p1a = points1(i, :);
    p1b = points1(i+1, :);

    for j = 1:n2-1
        p2a = points2(j, :);
        p2b = points2(j+1, :);

        % Direction vectors
        d1 = p1b - p1a;
        d2 = p2b - p2a;
        d12 = p1a - p2a;

        denom = d1(1)*d2(2) - d1(2)*d2(1);

        if abs(denom) < 1e-10
            continue;  % Parallel segments
        end

        t1 = (d2(1)*d12(2) - d2(2)*d12(1)) / denom;
        t2 = (d1(1)*d12(2) - d1(2)*d12(1)) / denom;

        if t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1
            pt = p1a + t1 * d1;
            return;
        end
    end
end
end

function lambda = intersect_ray_cylinder(ray, P_axis, axis_dir, R, lambda_ref)
% Find intersection of ray with cylinder
% Solves: ||cross(lambda*ray - P_axis, axis_dir)|| = R

ray = ray(:)';
P_axis = P_axis(:)';
axis_dir = axis_dir(:)';

u = cross(ray, axis_dir);
v = cross(P_axis, axis_dir);

A = dot(u, u);
B = -2 * dot(u, v);
C = dot(v, v) - R^2;

discriminant = B^2 - 4*A*C;

if discriminant < 0 || A < 1e-10
    lambda = lambda_ref;
    return;
end

sqrt_disc = sqrt(discriminant);
lambda1 = (-B + sqrt_disc) / (2*A);
lambda2 = (-B - sqrt_disc) / (2*A);

% Choose positive solution closest to reference
candidates = [lambda1, lambda2];
candidates = candidates(candidates > 0);

if isempty(candidates)
    lambda = lambda_ref;
else
    [~, idx] = min(abs(candidates - lambda_ref));
    lambda = candidates(idx);
end
end

% Submission.m - Integrated Architecture and Computer Vision Homework 2025-2026
% Project Title: 3D Geometric Reconstruction of a Cylindrical Vault
% Purpose: This script performs a full pipeline from image rectification to 3D recovery.

clear;
close all;
clc;

%% Part 0: Configuration and Setup
% In this section, we define the input paths for our image and feature set.
imageFile = 'San Maurizio.jpg';
featureFile = 'features_fixed.mat';

%% Part 1: Image and Feature Set Loading
% Here we load the primary image and the manually/automatically extracted features.
% These features include vertical lines, axis lines, and transversal lines which are
% essential for the vanishing point estimation.
fprintf('--- Loading Image and Feature Set ---\n');
img = imread(imageFile);
[rows, cols, ~] = size(img);
load(featureFile);
fprintf('Successfully loaded %d vertical, %d axis, and %d transversal lines.\n', ...
    length(lines_v), length(lines_axis), length(lines_trans));
fprintf('Identified %d arcs of type A and %d arcs of type B.\n', length(arcs_A), length(arcs_B));

%% Part 2: Geometric Computation of Vanishing Points
% The objective of this step is to find the vanishing points corresponding to the
% three principal orthogonal directions. We use a normalized coordinate system
% centered on the image to improve numerical stability during the SVD computation.
fprintf('\n--- Solving for Vanishing Points ---\n');

% We define the image center and a scaling factor for normalization purposes.
cx = cols / 2;
cy = rows / 2;
scale = max(cx, cy);

% Construction of the normalization matrix to map pixel space to normalized space.
T_norm_to_pixel = [scale, 0, cx; 0, scale, cy; 0, 0, 1];

% --- Directional Analysis: Vertical Vanishing Point ---
% We gather all vertical line segments and find their common intersection using SVD.
if length(lines_v) >= 2
    A_vert = zeros(length(lines_v), 3);
    for i = 1:length(lines_v)
        % Using our local function to compute the normalized line representation.
        line_hom = get_normalized_line(lines_v(i).p1, lines_v(i).p2, cx, cy, scale);
        A_vert(i, :) = line_hom';
    end
    [~, ~, V] = svd(A_vert);
    v_vert_normalized = V(:, end)';
    v_vert_normalized = v_vert_normalized / v_vert_normalized(3);
else
    error('Insufficient vertical lines (minimum 2 required).');
end

% --- Directional Analysis: Axis Vanishing Point ---
% This vanishing point represents the direction of the barrel vault's main axis.
% We consolidate the axis lines and the apical line to refine the estimation.
all_axis_lines = {};
for i = 1:length(lines_axis)
    all_axis_lines{end+1} = get_normalized_line(lines_axis(i).p1, lines_axis(i).p2, cx, cy, scale);
end
if ~isempty(line_apical)
    all_axis_lines{end+1} = get_normalized_line(line_apical.p1, line_apical.p2, cx, cy, scale);
end

A_axis = zeros(length(all_axis_lines), 3);
for i = 1:length(all_axis_lines)
    A_axis(i, :) = all_axis_lines{i}';
end
[~, ~, V] = svd(A_axis);
v_axis_normalized = V(:, end)';
v_axis_normalized = v_axis_normalized / v_axis_normalized(3);

% --- Directional Analysis: Transversal Vanishing Point ---
% This point corresponds to the direction perpendicular to both axis and vertical.
if length(lines_trans) >= 2
    A_trans = zeros(length(lines_trans), 3);
    for i = 1:length(lines_trans)
        line_hom = get_normalized_line(lines_trans(i).p1, lines_trans(i).p2, cx, cy, scale);
        A_trans(i, :) = line_hom';
    end
    [~, ~, V] = svd(A_trans);
    v_trans_normalized = V(:, end)';
    v_trans_normalized = v_trans_normalized / v_trans_normalized(3);
else
    error('Insufficient transversal lines (minimum 2 required).');
end

% We project the normalized vanishing points back into pixel coordinates for visualization.
v_vert_pixel = (T_norm_to_pixel * v_vert_normalized(:))';
v_vert_pixel = v_vert_pixel / v_vert_pixel(3);

v_axis_pixel = (T_norm_to_pixel * v_axis_normalized(:))';
v_axis_pixel = v_axis_pixel / v_axis_pixel(3);

v_trans_pixel = (T_norm_to_pixel * v_trans_normalized(:))';
v_trans_pixel = v_trans_pixel / v_trans_pixel(3);

% The vanishing line is computed as the line connecting two vanishing points.
% This line represents the horizon in the image plane for the specific orientation.
vanishing_line = cross(v_vert_pixel, v_trans_pixel);
vanishing_line = vanishing_line / norm(vanishing_line(1:2));

fprintf('Computed Vertical VP: (%.1f, %.1f)\n', v_vert_pixel(1), v_vert_pixel(2));
fprintf('Computed Axis VP:     (%.1f, %.1f)\n', v_axis_pixel(1), v_axis_pixel(2));
fprintf('Computed Trans VP:    (%.1f, %.1f)\n', v_trans_pixel(1), v_trans_pixel(2));

%% Part 3a: Projective and Affine Rectification
% This stage is dedicated to removing perspective distortion from the image.
% By identifying the vanishing line, we can apply a transformation that maps
% it back to infinity, effectively making parallel lines in the scene parallel
% in the rectified image.
fprintf('\n--- Part 3a: Executing Affine Rectification ---\n');

% Step 1: Affine rectification.
% We normalize the vanishing line and construct a homography that stabilizes the plane.
% Any line on the world plane maps to the vanishing line l = [l1, l2, l3] in the image.
% The transformation H_affine = [1 0 0; 0 1 0; l1 l2 l3] maps the vanishing line to [0 0 1].
l = vanishing_line(:) / vanishing_line(3);
H_affine = [1, 0, 0; 0, 1, 0; l(1), l(2), 1];

% Step 2: Basic Metric rectification (via Vanishing Points).
% In affine space, the vanishing line has been mapped to infinity.
% The vanishing points v_vert and v_trans now represent directions.
% Their transformed homogeneous coordinates are [x, y, 0]'.
v_vert_affine = H_affine * v_vert_pixel(:);
v_trans_affine = H_affine * v_trans_pixel(:);

% The direction vectors are extracted from the first two components.
u = v_vert_affine(1:2) / norm(v_vert_affine(1:2));
v = v_trans_affine(1:2) / norm(v_trans_affine(1:2));

% This is a simplified metric rectification that assumes the vertical and
% transversal directions define the new Cartesian frame.
M = [u, v];
S_metric_init = inv(M); % Maps u -> [1;0] and v -> [0;1]
H_metric_vp = [S_metric_init, [0; 0]; 0, 0, 1];

%% Part 3b: Metric Rectification via Independent Line Pairs (Automatic)
% This method recovers the metric structure of the plane by solving for the
% image of the dual conic C* using two pairs of orthogonal lines.
% We automatically try all permutations of v (vertical) and t (transversal)
% lines to find the best combination based on the condition number of S.
% Ref: Lecture G - Stratified Rectification from Orthogonal Lines.
fprintf('\n--- Part 3b: Metric Rectification (Automatic Permutation Search) ---\n');

% 1. VISUALIZE FEATURES
figure('Name', 'Part 3b: Automatic Orthogonal Line Pair Selection');
imshow(img); hold on;
title('Automatic Selection of Orthogonal Line Pairs (v-t permutations)');

% Plot features
plot_labeled_line_set(lines_v, 'b', 'v'); % Vertical in Blue
plot_labeled_line_set(lines_trans, 'r', 't'); % Transversal in Red
if ~isempty(lines_axis), plot_labeled_line_set(lines_axis, 'g', 'a'); end % Axis in Green
if ~isempty(line_apical), plot([line_apical.p1(1), line_apical.p2(1)], [line_apical.p1(2), line_apical.p2(2)], 'y-', 'LineWidth', 3); text(mean([line_apical.p1(1), line_apical.p2(1)]), mean([line_apical.p1(2), line_apical.p2(2)]), 'APIC', 'Color', 'y'); end

% 2. AUTOMATIC PERMUTATION SEARCH
% We need two pairs of orthogonal lines: (v_i, t_j) and (v_k, t_l)
% where i != k and j != l to ensure independence.
transform_line = @(L, H) cross((H*[L.p1, 1]')', (H*[L.p2, 1]')');

num_v = length(lines_v);
num_t = length(lines_trans);

fprintf('Testing all permutations of %d vertical x %d transversal lines...\n', num_v, num_t);

% Store results for all valid permutations
all_results = [];
result_idx = 0;

% Precompute vanishing point constraints for calibration (these don't change)
vps = {v_vert_pixel, v_axis_pixel, v_trans_pixel};
vp_pairs = [1 2; 1 3; 2 3];
A_vp = [];
for kk = 1:size(vp_pairs, 1)
    u_vp = vps{vp_pairs(kk, 1)};
    v_vp = vps{vp_pairs(kk, 2)};
    row_vp = [u_vp(1)*v_vp(1), u_vp(2)*v_vp(2), u_vp(1)*v_vp(3)+u_vp(3)*v_vp(1), u_vp(2)*v_vp(3)+u_vp(3)*v_vp(2), u_vp(3)*v_vp(3)];
    A_vp = [A_vp; row_vp];
end

% Generate all valid permutations: pick 2 distinct v lines and 2 distinct t lines
for v1_idx = 1:num_v
    for v2_idx = 1:num_v
        if v2_idx == v1_idx, continue; end % Must be different v lines

        for t1_idx = 1:num_t
            for t2_idx = 1:num_t
                if t2_idx == t1_idx, continue; end % Must be different t lines

                try
                    % Get line structures
                    L1 = lines_v(v1_idx);
                    M1 = lines_trans(t1_idx);
                    L2 = lines_v(v2_idx);
                    M2 = lines_trans(t2_idx);

                    % Transform to affine space
                    l1_a = transform_line(L1, H_affine);
                    m1_a = transform_line(M1, H_affine);
                    l2_a = transform_line(L2, H_affine);
                    m2_a = transform_line(M2, H_affine);

                    % Normalize lines
                    l1_a = l1_a / l1_a(3); m1_a = m1_a / m1_a(3);
                    l2_a = l2_a / l2_a(3); m2_a = m2_a / m2_a(3);

                    % Build constraint matrix for S
                    A_metric = [l1_a(1)*m1_a(1), (l1_a(1)*m1_a(2) + l1_a(2)*m1_a(1)), l1_a(2)*m1_a(2);
                        l2_a(1)*m2_a(1), (l2_a(1)*m2_a(2) + l2_a(2)*m2_a(1)), l2_a(2)*m2_a(2)];

                    % Check rank
                    if rank(A_metric) < 2, continue; end

                    % Solve for S parameters
                    [~, ~, V_s] = svd(A_metric);
                    s_params = V_s(:, end);
                    S = [s_params(1), s_params(2); s_params(2), s_params(3)];

                    % Ensure positive definiteness
                    if det(S) < 0 || S(1,1) < 0, S = -S; end
                    if det(S) <= 0 || S(1,1) <= 0, continue; end

                    % Try Cholesky decomposition
                    L_chol = chol(S, 'lower');
                    H_metric_test = [inv(L_chol), [0; 0]; 0, 0, 1];

                    % Compute quality metrics
                    cond_S = cond(S);
                    metric_aspect_ratio = L_chol(2,2) / L_chol(1,1);

                    % === COMPUTE CAMERA CALIBRATION FOR THIS PERMUTATION ===
                    H_rectify_test = H_metric_test * H_affine;
                    H_world_to_img = inv(H_rectify_test);
                    h1 = H_world_to_img(:, 1);
                    h2 = H_world_to_img(:, 2);

                    % Build IAC constraints: VP orthogonality + rectification constraints
                    A_iac_test = A_vp;

                    % Constraint: Rectified axes orthogonal (h1' * omega * h2 = 0)
                    row_ortho = [h1(1)*h2(1), h1(2)*h2(2), h1(1)*h2(3)+h1(3)*h2(1), h1(2)*h2(3)+h1(3)*h2(2), h1(3)*h2(3)];
                    A_iac_test = [A_iac_test; row_ortho];

                    % Constraint: Equal scale (h1'*omega*h1 = h2'*omega*h2)
                    t1_iac = [h1(1)*h1(1), h1(2)*h1(2), 2*h1(1)*h1(3), 2*h1(2)*h1(3), h1(3)*h1(3)];
                    t2_iac = [h2(1)*h2(1), h2(2)*h2(2), 2*h2(1)*h2(3), 2*h2(2)*h2(3), h2(3)*h2(3)];
                    A_iac_test = [A_iac_test; (t1_iac - t2_iac)];

                    % Solve for omega
                    [~, ~, V_iac] = svd(A_iac_test);
                    x_iac = V_iac(:, end);
                    omega_test = [x_iac(1), 0, x_iac(3); 0, x_iac(2), x_iac(4); x_iac(3), x_iac(4), x_iac(5)];

                    % Force positive definiteness
                    if det(omega_test) < 0, omega_test = -omega_test; end

                    % Extract K via Cholesky
                    C_iac = chol(omega_test, 'upper');
                    K_test = inv(C_iac);
                    K_test = K_test / K_test(3,3);
                    if K_test(1,1) < 0, K_test = -K_test; end

                    % Extract calibration parameters
                    fx = K_test(1,1);
                    fy = K_test(2,2);
                    cam_aspect_ratio = fx / fy;
                    pp_x = K_test(1,3);
                    pp_y = K_test(2,3);

                    % Store result
                    result_idx = result_idx + 1;
                    all_results(result_idx).v1 = v1_idx;
                    all_results(result_idx).t1 = t1_idx;
                    all_results(result_idx).v2 = v2_idx;
                    all_results(result_idx).t2 = t2_idx;
                    all_results(result_idx).cond_S = cond_S;
                    all_results(result_idx).aspect_ratio = metric_aspect_ratio;
                    all_results(result_idx).S = S;
                    all_results(result_idx).H_metric = H_metric_test;
                    all_results(result_idx).K = K_test;
                    all_results(result_idx).fx = fx;
                    all_results(result_idx).fy = fy;
                    all_results(result_idx).cam_ar = cam_aspect_ratio;
                    all_results(result_idx).pp_x = pp_x;
                    all_results(result_idx).pp_y = pp_y;

                catch
                    % Skip invalid combinations
                    continue;
                end
            end
        end
    end
end

fprintf('\nFound %d valid permutations.\n', length(all_results));

% 3. DISPLAY ALL RESULTS COMPARISON
if ~isempty(all_results)
    % Compute statistics across all permutations
    cond_values = [all_results.cond_S];
    aspect_values = [all_results.aspect_ratio];
    cam_ar_values = [all_results.cam_ar];
    fx_values = [all_results.fx];
    fy_values = [all_results.fy];
    pp_x_values = [all_results.pp_x];
    pp_y_values = [all_results.pp_y];

    fprintf('\n========== PERMUTATION STATISTICS ==========\n');
    fprintf('Total valid permutations: %d\n', length(all_results));
    fprintf('\nCondition Number Statistics:\n');
    fprintf('  Min:    %.4f\n', min(cond_values));
    fprintf('  Max:    %.4f\n', max(cond_values));
    fprintf('  Mean:   %.4f\n', mean(cond_values));
    fprintf('  Median: %.4f\n', median(cond_values));
    fprintf('  Std:    %.4f\n', std(cond_values));

    fprintf('\nMetric Aspect Ratio (from S) Statistics:\n');
    fprintf('  Min:    %.4f\n', min(aspect_values));
    fprintf('  Max:    %.4f\n', max(aspect_values));
    fprintf('  Mean:   %.4f\n', mean(aspect_values));
    fprintf('  Median: %.4f\n', median(aspect_values));
    fprintf('  Std:    %.4f\n', std(aspect_values));

    fprintf('\n========== CAMERA CALIBRATION STATISTICS ==========\n');
    fprintf('\nCamera Aspect Ratio (fx/fy) Statistics:\n');
    fprintf('  Min:    %.4f\n', min(cam_ar_values));
    fprintf('  Max:    %.4f\n', max(cam_ar_values));
    fprintf('  Mean:   %.4f\n', mean(cam_ar_values));
    fprintf('  Median: %.4f\n', median(cam_ar_values));
    fprintf('  Std:    %.4f\n', std(cam_ar_values));

    fprintf('\nFocal Length fx Statistics:\n');
    fprintf('  Min:    %.2f\n', min(fx_values));
    fprintf('  Max:    %.2f\n', max(fx_values));
    fprintf('  Mean:   %.2f\n', mean(fx_values));
    fprintf('  Median: %.2f\n', median(fx_values));
    fprintf('  Std:    %.2f\n', std(fx_values));

    fprintf('\nFocal Length fy Statistics:\n');
    fprintf('  Min:    %.2f\n', min(fy_values));
    fprintf('  Max:    %.2f\n', max(fy_values));
    fprintf('  Mean:   %.2f\n', mean(fy_values));
    fprintf('  Median: %.2f\n', median(fy_values));
    fprintf('  Std:    %.2f\n', std(fy_values));

    fprintf('\nPrincipal Point X Statistics:\n');
    fprintf('  Min:    %.2f\n', min(pp_x_values));
    fprintf('  Max:    %.2f\n', max(pp_x_values));
    fprintf('  Mean:   %.2f\n', mean(pp_x_values));
    fprintf('  Median: %.2f\n', median(pp_x_values));
    fprintf('  Std:    %.2f\n', std(pp_x_values));

    fprintf('\nPrincipal Point Y Statistics:\n');
    fprintf('  Min:    %.2f\n', min(pp_y_values));
    fprintf('  Max:    %.2f\n', max(pp_y_values));
    fprintf('  Mean:   %.2f\n', mean(pp_y_values));
    fprintf('  Median: %.2f\n', median(pp_y_values));
    fprintf('  Std:    %.2f\n', std(pp_y_values));

    % Sort by condition number (lower is better)
    [~, sort_idx] = sort(cond_values);

    % ========== TOP 200 CAMERA CALIBRATION RESULTS ==========
    num_to_show = min(200, length(all_results));
    fprintf('\n========== TOP %d PERMUTATIONS: CAMERA CALIBRATION ==========\n', num_to_show);
    fprintf('%-5s | %-14s | %-8s | %-8s | %-8s | %-10s | %-10s\n', ...
        'Rank', 'Lines', 'fx', 'fy', 'AR(fx/fy)', 'PP_x', 'PP_y');
    fprintf('%s\n', repmat('-', 1, 85));

    for i = 1:num_to_show
        idx = sort_idx(i);
        r = all_results(idx);
        line_str = sprintf('v%d-t%d,v%d-t%d', r.v1, r.t1, r.v2, r.t2);
        fprintf('%-5d | %-14s | %-8.1f | %-8.1f | %-8.4f | %-10.2f | %-10.2f\n', ...
            i, line_str, r.fx, r.fy, r.cam_ar, r.pp_x, r.pp_y);
    end

    % Show distribution of quality categories
    excellent_count = sum(cond_values < 10 & abs(aspect_values - 1) < 0.3);
    good_count = sum(cond_values < 50 & abs(aspect_values - 1) < 0.5) - excellent_count;
    acceptable_count = sum(cond_values < 100) - excellent_count - good_count;
    poor_count = length(all_results) - excellent_count - good_count - acceptable_count;

    fprintf('\n========== QUALITY DISTRIBUTION ==========\n');
    fprintf('Excellent:   %3d (%5.1f%%)\n', excellent_count, 100*excellent_count/length(all_results));
    fprintf('Good:        %3d (%5.1f%%)\n', good_count, 100*good_count/length(all_results));
    fprintf('Acceptable:  %3d (%5.1f%%)\n', acceptable_count, 100*acceptable_count/length(all_results));
    fprintf('Poor:        %3d (%5.1f%%)\n', poor_count, 100*poor_count/length(all_results));

    % Show BEST vs WORST comparison
    fprintf('\n========== BEST vs WORST COMPARISON ==========\n');
    best_r = all_results(sort_idx(1));
    worst_r = all_results(sort_idx(end));

    fprintf('BEST:  v%d-t%d, v%d-t%d\n', best_r.v1, best_r.t1, best_r.v2, best_r.t2);
    fprintf('       Cond=%.4f | AR=%.4f | fx=%.1f | fy=%.1f | PP=(%.1f, %.1f)\n', ...
        best_r.cond_S, best_r.aspect_ratio, best_r.fx, best_r.fy, best_r.pp_x, best_r.pp_y);
    fprintf('WORST: v%d-t%d, v%d-t%d\n', worst_r.v1, worst_r.t1, worst_r.v2, worst_r.t2);
    fprintf('       Cond=%.4f | AR=%.4f | fx=%.1f | fy=%.1f | PP=(%.1f, %.1f)\n', ...
        worst_r.cond_S, worst_r.aspect_ratio, worst_r.fx, worst_r.fy, worst_r.pp_x, worst_r.pp_y);
    fprintf('Ratio (worst/best cond): %.2fx\n', worst_r.cond_S / best_r.cond_S);

    % 4. SELECT BEST RESULT
    best_idx = sort_idx(1);
    best_result = all_results(best_idx);
    H_metric = best_result.H_metric;

    fprintf('\n========== SELECTED BEST PERMUTATION ==========\n');
    fprintf('Pair 1: v%d + t%d\n', best_result.v1, best_result.t1);
    fprintf('Pair 2: v%d + t%d\n', best_result.v2, best_result.t2);
    fprintf('Condition Number: %.4f\n', best_result.cond_S);
    fprintf('Aspect Ratio: %.4f\n', best_result.aspect_ratio);
    fprintf('Metric rectification successful (automatic selection).\n');

else
    % Fallback to VP-based method if no valid permutation found
    warning('No valid v-t permutation found. Falling back to VP-based metric rectification.');
    H_metric = H_metric_vp;
end

% Final combined homography
H_rectify = H_metric * H_affine;

% --- Visualization and Output Generation ---
% In order to visualize the result properly, we need to compute the output bounds.
[xx, yy] = meshgrid(linspace(1, cols, 40), linspace(1, rows, 40));
grid_points = [xx(:), yy(:), ones(numel(xx), 1)]';

grid_vals = H_affine(3, :) * grid_points;
side = sign(mean(grid_vals));
if side == 0, side = 1; end

% We apply a validity mask to focus on the points that are not near the vanishing line.
margin = 0.15 * abs(max(grid_vals) - min(grid_vals));
valid_mask = (sign(grid_vals) == side) & (abs(grid_vals) > margin);
if sum(valid_mask) < 10, valid_mask = (sign(grid_vals) == side); end

pts_transformed = H_rectify * grid_points(:, valid_mask);
pts_transformed = pts_transformed(1:2, :) ./ pts_transformed(3, :);

% Bounding box estimation for the rectified image.
min_x = min(pts_transformed(1, :)); max_x = max(pts_transformed(1, :));
min_y = min(pts_transformed(2, :)); max_y = max(pts_transformed(2, :));

rect_width = max_x - min_x; rect_height = max_y - min_y;

% Scaling the rectified result to a standard presentation width (2000 pixels).
target_width = 2000;
scale_factor = target_width / rect_width;

% Final transformation including translation to positive coordinates.
T_output = [scale_factor, 0, -scale_factor*min_x + 1; ...
    0, scale_factor, -scale_factor*min_y + 1; ...
    0, 0, 1];

H_final = T_output * H_rectify;

output_size = [ceil(scale_factor * rect_height), target_width];
img_rectified = imwarp(img, projective2d(H_final'), 'OutputView', imref2d(output_size));

fprintf('Rectified image generated with size: %d x %d\n', output_size(2), output_size(1));

%% Part 4: Camera Intrinsic Calibration (Robust IAC Method)
% We estimate K without assuming fx=fy or centered principal point.
% We use the Image of the Absolute Conic (omega).
% Assumption: Zero skew (given), but unknown Aspect Ratio and Principal Point.

fprintf('\n--- Computing Camera Calibration matrix K (General Method) ---\n');

% 1. Collect Constraints
% We define omega as symmetric with zero skew:
% omega = [x1  0  x3]
%         [ 0 x2  x4]
%         [x3 x4  x5]
% The vector of unknowns is x = [x1, x2, x3, x4, x5]'.

A_iac = [];

% Constraint Set A: Vanishing Points Orthogonality
% vi' * omega * vj = 0 for orthogonal directions
vps = {v_vert_pixel, v_axis_pixel, v_trans_pixel};
pairs = [1 2; 1 3; 2 3]; % (Vert-Axis), (Vert-Trans), (Axis-Trans)

for k = 1:size(pairs, 1)
    u = vps{pairs(k, 1)};
    v = vps{pairs(k, 2)};

    % Expansion of u' * omega * v = 0 with zero skew structure
    % x1(u1v1) + x2(u2v2) + x3(u1v3+u3v1) + x4(u2v3+u3v2) + x5(u3v3) = 0
    row = [u(1)*v(1), ...
        u(2)*v(2), ...
        u(1)*v(3) + u(3)*v(1), ...
        u(2)*v(3) + u(3)*v(2), ...
        u(3)*v(3)];
    A_iac = [A_iac; row];
end

% Constraint Set B: Scene Geometry from Rectification
% We use the rectification homography H_rectify computed in Part 3.
% This maps Image -> World (Metric).
% Therefore H_inv = inv(H_rectify) maps World -> Image.
% The columns of H_inv represent the World X and Y axes in the Image.
H_img_to_world = H_rectify; % From your Part 3
H_world_to_img = inv(H_img_to_world);

h1 = H_world_to_img(:, 1); % Image of World X axis
h2 = H_world_to_img(:, 2); % Image of World Y axis

% Constraint 4: Rectified axes are orthogonal in 3D (h1' * omega * h2 = 0)
u = h1; v = h2;
row_ortho = [u(1)*v(1), ...
    u(2)*v(2), ...
    u(1)*v(3) + u(3)*v(1), ...
    u(2)*v(3) + u(3)*v(2), ...
    u(3)*v(3)];
A_iac = [A_iac; row_ortho];

% Constraint 5: Rectified axes have equal scale (h1' * omega * h1 = h2' * omega * h2)
% Equivalent to: h1' * omega * h1 - h2' * omega * h2 = 0
% Term 1 (h1, h1)
t1 = [h1(1)*h1(1), h1(2)*h1(2), 2*h1(1)*h1(3), 2*h1(2)*h1(3), h1(3)*h1(3)];
% Term 2 (h2, h2)
t2 = [h2(1)*h2(1), h2(2)*h2(2), 2*h2(1)*h2(3), 2*h2(2)*h2(3), h2(3)*h2(3)];

A_iac = [A_iac; (t1 - t2)];

% 2. Solve for omega using SVD
[~, ~, V_calib] = svd(A_iac);
x = V_calib(:, end);

% Reconstruct omega matrix
omega = [x(1), 0,    x(3); ...
    0,    x(2), x(4); ...
    x(3), x(4), x(5)];

% 3. Extract K using Cholesky Factorization
% omega = inv(K * K')
% K_inv = cholesky(omega)
try
    % Force positive definiteness if flip occurred during SVD
    if det(omega) < 0
        omega = -omega;
    end

    C = chol(omega, 'upper'); % C such that C'*C = omega
    K = inv(C);               % K is the inverse of the Cholesky factor

    % Normalize K so K(3,3) = 1
    K = K / K(3,3);

    % Ensure positive focal lengths
    if K(1,1) < 0, K = -K; end

    fprintf('Robust Calibration Successful.\n');
    fprintf('f_x = %.2f\n', K(1,1));
    fprintf('f_y = %.2f\n', K(2,2));
    fprintf('Aspect Ratio = %.4f\n', K(1,1)/K(2,2));
    fprintf('Principal Point = (%.2f, %.2f)\n', K(1,3), K(2,3));

catch
    warning('Cholesky decomposition failed. Matrix omega may not be positive definite due to feature noise.');
    % Fallback to simplified model ONLY if math fails significantly
    fprintf('Falling back to simplified calibration for stability.\n');
    K = [2000, 0, cx; 0, 2000, cy; 0, 0, 1];
end


%% Part 5: Determination of Nodal Points
% Nodal points are the critical intersection points of the vault's arcs.
% Finding these points allows us to define the skeleton of the reconstruction.
fprintf('\n--- Searching for Nodal Intersection Points ---\n');

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

% Apply known index shift correction (B_index offset from A_index).
shift = 1;
nodal_points(:, 4) = nodal_points(:, 4) - shift;

% Print apical nodes (where A_index == B_index).
fprintf('\n--- Apical Nodal Points (on cylinder apex) ---\n');
for i = 1:size(nodal_points, 1)
    if nodal_points(i, 3) == nodal_points(i, 4)
        fprintf('APICAL Node N(%d,%d) at (%.0f, %.0f)\n', ...
            nodal_points(i, 3), nodal_points(i, 4), nodal_points(i, 1), nodal_points(i, 2));
    end
end

% Select a specific non-apical nodal point N(1,2) for demonstration.
fprintf('\n--- Non-Apical Nodal Point (on cylinder surface) ---\n');
non_apical_idx = find(nodal_points(:,3) == 1 & nodal_points(:,4) == 2, 1);
np = nodal_points(non_apical_idx, :);
fprintf('Non-apical Node N(%d,%d) at image coords: (%.1f, %.1f)\n', np(3), np(4), np(1), np(2));
np_rect = H_final * [np(1:2), 1]';
np_rect = np_rect(1:2) / np_rect(3);
fprintf('Same point in rectified plane: (%.1f, %.1f)\n', np_rect(1), np_rect(2));


%% Part 6: Geometric 3D Reconstruction and Radius Estimation
% The final computational stage involves lifting our 2D measurements into
% 3D space. We use the calibrated K matrix and the vanishing point directions
% to define rays in space.
fprintf('\n--- Initiating 3D Geometric Reconstruction ---\n');

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

% Compute reference depth using two apical nodes via linear triangulation.
% The key insight is that P2 = P1 + delta_d * axis_direction.
% Substituting P1 = lambda1*ray1 and P2 = lambda2*ray2:
%   lambda2 * ray2 = lambda1 * ray1 + delta_d * axis_direction
% Rearranging: [ray1, -ray2] * [lambda1; lambda2] = -delta_d * axis_direction
n1_2d = apical_nodes(1, 1:2);
n2_2d = apical_nodes(2, 1:2);

ray1 = K_inv * [n1_2d, 1]';
ray1 = ray1 / norm(ray1);

ray2 = K_inv * [n2_2d, 1]';
ray2 = ray2 / norm(ray2);

idx1 = apical_nodes(1, 3);
idx2 = apical_nodes(2, 3);
delta_d = (idx2 - idx1) * d;

% Construct the 3x2 linear system and solve using least squares.
A_tri = [ray1(:), -ray2(:)];
b_tri = -delta_d * axis_direction(:);
lambdas_init = A_tri \ b_tri;
lambda1 = lambdas_init(1);
lambda2 = lambdas_init(2);

fprintf('Triangulated depths: lambda1 = %.3f, lambda2 = %.3f\n', lambda1, lambda2);

% Use lambda1 as the reference depth for subsequent calculations.
lambda_ref = lambda1;
a1 = dot(ray1, axis_direction);  % Needed for non-apical node depth scaling
fprintf('Reference depth (lambda_ref): %.3f\n', lambda_ref);

% 3D position of the first apical node (apex reference).
P_apex = lambda1 * ray1(:)';


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

% We estimate the cylinder radius by analyzing the distance of nodal points
% from the apex reference point. Outliers are removed via median filtering.
R_cylinder = median(radius_estimates);
fprintf('Estimated Cylinder Radius R: %.3f\n', R_cylinder);

% Cylinder axis localization (positioning the axis below the apex).
P_axis = P_apex - R_cylinder * vert_direction(:)';

% Documentation of the cylinder axis localization relative to the camera.
fprintf('\n--- Cylinder Axis Spatial Localization ---\n');
fprintf('The Axis passes through: [%.4f, %.4f, %.4f]\n', P_axis);
fprintf('Axis Direction Vector:    [%.4f, %.4f, %.4f]\n', axis_direction);

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

% Verify apical node distance using triangulated depths.
P1_apical = lambda1 * ray1(:)';
P2_apical = lambda2 * ray2(:)';
apical_distance = norm(P2_apical - P1_apical);
apical_axis_dist = abs(dot(P2_apical - P1_apical, axis_direction));
fprintf('Apical node 3D distance: %.3f (along axis: %.3f, target: %.1f)\n', apical_distance, apical_axis_dist, abs(delta_d));


% Print 3D coordinates of a dozen points from one arc (requirement 5)
fprintf('\n--- 3D Coordinates of 12 Points from Arc A1 (diagonal arc) ---\n');
arc1_pts = points_3D{1}.pts;
num_to_show = min(12, size(arc1_pts, 1));
fprintf('Point |     X     |     Y     |     Z     \n');
fprintf('------+-----------+-----------+-----------\n');
for i = 1:num_to_show
    fprintf('  %2d  | %9.4f | %9.4f | %9.4f\n', i, arc1_pts(i, 1), arc1_pts(i, 2), arc1_pts(i, 3));
end

%% Part 7: Data Visualization and Result Analysis
% In this final part, we generate several figures to visualize the
% reconstruction quality and the geometric relationships between
% vanishing points, nodal points, and the 3D arcs.
fprintf('\n--- Generating Comprehensive Visualizations ---\n');

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

%% Part 8: Termination and Cleanup
fprintf('\n--- Completed Successfully ---\n');

%% Helper Functions and Auxiliary Routines
% These functions encapsulate specific geometric operations used throughout the script.

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

function line_hom = get_normalized_line(p1, p2, cx, cy, scale)
% Computes a homogeneous line representation in normalized coordinates.
% This improves numerical stability for vanishing point estimation.
p1_norm = [(p1(1)-cx)/scale, (p1(2)-cy)/scale, 1];
p2_norm = [(p2(1)-cx)/scale, (p2(2)-cy)/scale, 1];
line_hom = cross(p1_norm, p2_norm)';
end

function out = ternary(cond, a, b)
% Simple ternary operator implementation
if cond, out = a; else, out = b; end
end

function plot_labeled_line_set(lines, color, prefix)
% Helper function to plot and label a set of lines
for i = 1:length(lines)
    plot([lines(i).p1(1), lines(i).p2(1)], [lines(i).p1(2), lines(i).p2(2)], color, 'LineWidth', 2);
    text(mean([lines(i).p1(1), lines(i).p2(1)]), mean([lines(i).p1(2), lines(i).p2(2)]), ...
        sprintf('%s%d', prefix, i), 'Color', color, 'FontSize', 10, 'FontWeight', 'bold');
end
end

function line = get_line_by_idx(idx, lines_v, lines_trans, lines_axis, line_apical)
% Local function to map a string index (e.g. 'v1', 't6') to the correct line structure.
idx = lower(idx);
num = str2double(idx(2:end));

if startsWith(idx, 'v') || startsWith(idx, 'l')
    line = lines_v(num);
elseif startsWith(idx, 't')
    line = lines_trans(num);
elseif startsWith(idx, 'a')
    line = lines_axis(num);
else
    line = line_apical;
end
end

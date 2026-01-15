% IACV Homework 2025-2026
% 3D Reconstruction of a Cylindric Vault
% Author: [Your Name]
% Date: 2026-01-03

clear; close all; clc;

%% Setup and Parameters
imageFile = 'San Maurizio.jpg';
featureFile = 'features_fixed copy.mat';

%% 1. Load Image
if exist(imageFile, 'file')
    img = imread(imageFile);
else
    error('Image file %s not found. Please make sure it is in the current directory.', imageFile);
end

[rows, cols, ~] = size(img);

%% 2. Feature Extraction (Smart Loading Version)
if ~exist('lines_v', 'var'), lines_v = []; end
if ~exist('lines_axis', 'var'), lines_axis = []; end
if ~exist('line_apical', 'var'), line_apical = []; end
if ~exist('lines_trans', 'var'), lines_trans = []; end
if ~exist('arcs_A', 'var'), arcs_A = {}; end
if ~exist('arcs_B', 'var'), arcs_B = {}; end

needs_save = false;

% Load primary features (lines, etc) from "copy" file
if exist('features_fixed copy.mat', 'file')
    fprintf('Loading primary features from features_fixed copy.mat...\n');
    load('features_fixed copy.mat');
end

% Overwrite specifically the arc features from "fixed" file as requested
if exist('features_fixed.mat', 'file')
    fprintf('Loading arc points from features_fixed.mat...\n');
    arc_data = load('features_fixed.mat', 'arcs_A', 'arcs_B');
    if isfield(arc_data, 'arcs_A'), arcs_A = arc_data.arcs_A; end
    if isfield(arc_data, 'arcs_B'), arcs_B = arc_data.arcs_B; end
end

repickList = [];
if exist(featureFile, 'file')
    prompt = {['Select features to RE-PICK (type numbers e.g. "1 3 5" or leave empty to KEEP ALL):', newline, ...
        '1. Vertical (Black)', newline, ...
        '2. Axis-Parallel (Green)', newline, ...
        '3. Apical (Yellow)', newline, ...
        '4. Transversal (White)', newline, ...
        '5. Arcs Family A (Cyan)', newline, ...
        '6. Arcs Family B (Magenta)']};
    answer = inputdlg(prompt, 'Refine Features', [1 100]);

    if ~isempty(answer) && ~isempty(answer{1})
        repickList = str2num(answer{1});
    end
end

if isempty(lines_v) || ismember(1, repickList)
    fprintf('1. Select VERTICAL lines (Black).\n');
    lines_v = select_lines(img, 'Select Vertical Lines (Black). Enter to finish.', 'k');
    needs_save = true;
end

if isempty(lines_axis) || ismember(2, repickList)
    fprintf('2. Select AXIS-PARALLEL lines (Green).\n');
    lines_axis = select_lines(img, 'Select Axis-Parallel Lines (Green). Enter to finish.', 'g');
    needs_save = true;
end

if isempty(line_apical) || ismember(3, repickList)
    fprintf('3. Select THE APICAL line (Yellow).\n');
    temp_apical = select_lines(img, 'Select THE Apical Line (Yellow). Select ONE and Enter.', 'y');
    if ~isempty(temp_apical), line_apical = temp_apical(1); end
    needs_save = true;
end

if isempty(lines_trans) || ismember(4, repickList)
    fprintf('4. Select TRANSVERSAL lines (White).\n');
    lines_trans = select_lines(img, 'Select Transversal Lines (White). Enter to finish.', 'w');
    needs_save = true;
end

if isempty(arcs_A) || ismember(5, repickList)
    fprintf('5. Select points for DIAGONAL ARCS (Family A - Cyan).\n');
    arcs_A = select_arcs(img, 'Select Family A Arcs (Up-Right /). Enter on empty to Finish');
    needs_save = true;
end

if isempty(arcs_B) || ismember(6, repickList)
    fprintf('6. Select points for DIAGONAL ARCS (Family B - Magenta).\n');
    arcs_B = select_arcs(img, 'Select Family B Arcs (Up-Left \). Enter twice to Finish.');
    needs_save = true;
end

if needs_save
    save(featureFile, 'lines_v', 'lines_axis', 'line_apical', 'lines_trans', 'arcs_A', 'arcs_B');
    fprintf('Features updated and saved to %s\n', featureFile);
end

figure(1); imshow(img); title('All Features Loaded'); hold on;
if ~isempty(lines_v), for i=1:length(lines_v), plot([lines_v(i).p1(1) lines_v(i).p2(1)], [lines_v(i).p1(2) lines_v(i).p2(2)], 'k-', 'LineWidth', 2); end; end
if ~isempty(lines_axis), for i=1:length(lines_axis), plot([lines_axis(i).p1(1) lines_axis(i).p2(1)], [lines_axis(i).p1(2) lines_axis(i).p2(2)], 'g-', 'LineWidth', 2); end; end
if ~isempty(line_apical), plot([line_apical.p1(1) line_apical.p2(1)], [line_apical.p1(2) line_apical.p2(2)], 'y-', 'LineWidth', 3); end
if ~isempty(lines_trans), for i=1:length(lines_trans), plot([lines_trans(i).p1(1) lines_trans(i).p2(1)], [lines_trans(i).p1(2) lines_trans(i).p2(2)], 'w-', 'LineWidth', 2); end; end
for i=1:length(arcs_A), pts=arcs_A{i}; plot(pts(:,1), pts(:,2), 'c.-', 'MarkerSize', 10); end
for i=1:length(arcs_B), pts=arcs_B{i}; plot(pts(:,1), pts(:,2), 'm.-', 'MarkerSize', 10); end

%% 3. Vanishing Points and Vanishing Line
fprintf('Computing Vanishing Points (NORMALIZED)...\n');

cx = cols / 2;
cy = rows / 2;
fprintf('  Normalizing data: shift to center (%.1f, %.1f)\n', cx, cy);

scale = max(cx, cy);
fprintf('  Scaling factor: %.2f\n', scale);

T_n2p = [scale 0 cx; 0 scale cy; 0 0 1];

get_line_n = @(p1, p2) cross([(p1(1)-cx)/scale, (p1(2)-cy)/scale, 1], ...
    [(p2(1)-cx)/scale, (p2(2)-cy)/scale, 1])';

v_vert_n = [];
if ~isempty(lines_v)
    lines_hom = {};
    for i=1:length(lines_v), lines_hom{i} = get_line_n(lines_v(i).p1, lines_v(i).p2); end
    v_vert_n = compute_vanishing_point(lines_hom);
    fprintf('v_vert_n = [%.4f, %.4f, %.4f]\n', v_vert_n);
end

v_axis_n = [];
lines_hom = {};
for i=1:length(lines_axis), lines_hom{end+1} = get_line_n(lines_axis(i).p1, lines_axis(i).p2); end
if ~isempty(line_apical), lines_hom{end+1} = get_line_n(line_apical.p1, line_apical.p2); end

if length(lines_hom) >= 2
    v_axis_n = compute_vanishing_point(lines_hom);
    fprintf('v_axis_n = [%.4f, %.4f, %.4f]\n', v_axis_n);
end

v_trans_n = [];
if ~isempty(lines_trans)
    lines_hom = {};
    for i=1:length(lines_trans), lines_hom{i} = get_line_n(lines_trans(i).p1, lines_trans(i).p2); end
    v_trans_n = compute_vanishing_point(lines_hom);
    fprintf('v_trans_n = [%.4f, %.4f, %.4f]\n', v_trans_n);
end

if ~isempty(v_vert_n), v_vert = (T_n2p * v_vert_n')'; else, v_vert = []; end
if ~isempty(v_axis_n), v_axis = (T_n2p * v_axis_n')'; else, v_axis = []; end
if ~isempty(v_trans_n), v_trans = (T_n2p * v_trans_n')'; else, v_trans = []; end

fprintf('Denormalized VPs (Pixel Coords):\n');
if ~isempty(v_vert), fprintf('  v_vert  = [%.4f, %.4f, %.4f]\n', v_vert); end
if ~isempty(v_axis), fprintf('  v_axis  = [%.4f, %.4f, %.4f]\n', v_axis); end
if ~isempty(v_trans), fprintf('  v_trans = [%.4f, %.4f, %.4f]\n', v_trans); end

% Compute the vanishing line of the plane perpendicular to the cylinder axis.
% This plane contains the vertical direction and the transversal direction.
if ~isempty(v_vert) && ~isempty(v_trans)
    l_inf_perp = cross(v_vert, v_trans);
    l_inf_perp = l_inf_perp / norm(l_inf_perp(1:2));
    fprintf('Vanishing line l_inf_perp = [%.4f, %.4f, %.4f]\n', l_inf_perp);
else
    l_inf_perp = [];
    warning('Missing vertical or transversal vanishing point. Cannot compute vanishing line.');
end

if ~isempty(v_vert) || ~isempty(v_axis) || ~isempty(v_trans)
    figure(2); clf; imshow(img); title('Vanishing Points & Line Extensions'); hold on;
    all_pts = [1, 1; cols, rows];

    % Plot vanishing points
    if ~isempty(v_vert)
        plot(v_vert(1)/v_vert(3), v_vert(2)/v_vert(3), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
        text(v_vert(1)/v_vert(3), v_vert(2)/v_vert(3), '  v_{vert}', 'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
    end
    if ~isempty(v_axis)
        plot(v_axis(1)/v_axis(3), v_axis(2)/v_axis(3), 'go', 'MarkerSize', 12, 'LineWidth', 2);
        text(v_axis(1)/v_axis(3), v_axis(2)/v_axis(3), '  v_{axis}', 'Color', 'g', 'FontSize', 10, 'FontWeight', 'bold');
    end
    if ~isempty(v_trans)
        plot(v_trans(1)/v_trans(3), v_trans(2)/v_trans(3), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
        text(v_trans(1)/v_trans(3), v_trans(2)/v_trans(3), '  v_{trans}', 'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
    end

    % Extensions to the vanishing points
    if ~isempty(v_vert)
        vp = v_vert(1:2)/v_vert(3); for i=1:length(lines_v), plot([lines_v(i).p1(1), vp(1)], [lines_v(i).p1(2), vp(2)], 'k:', 'LineWidth', 0.5); end
    end
    if ~isempty(v_axis)
        vp = v_axis(1:2)/v_axis(3); for i=1:length(lines_axis), plot([lines_axis(i).p1(1), vp(1)], [lines_axis(i).p1(2), vp(2)], 'g:', 'LineWidth', 0.5); end
    end
    if ~isempty(v_trans)
        vp = v_trans(1:2)/v_trans(3); for i=1:length(lines_trans), plot([lines_trans(i).p1(1), vp(1)], [lines_trans(i).p1(2), vp(2)], 'r:', 'LineWidth', 0.5); end
    end

    % Draw the vanishing line
    if ~isempty(l_inf_perp)
        a = l_inf_perp(1); b = l_inf_perp(2); c = l_inf_perp(3);

        % Calculate limits for drawing
        xl = [1, cols];
        yl = [1, rows];

        % Use a very large bounding box for drawing the line
        draw_x = [-cols*5, cols*5];
        if abs(b) > 1e-9
            draw_y = -(a*draw_x + c)/b;
            plot(draw_x, draw_y, 'b-', 'LineWidth', 2.5);
            text(draw_x(2), draw_y(2), '  l_{\infty\perp}', 'Color', 'b', 'FontSize', 14, 'FontWeight', 'bold');
        else
            x_v = -c/a;
            plot([x_v, x_v], [-rows*5, rows*5], 'b-', 'LineWidth', 2.5);
        end
    end

    % Adjust view
    set(gca, 'Visible', 'on', 'Layer', 'top');
    grid on; grid minor;
    set(gca, 'Color', [0.9 0.9 0.9]);
end

%% 4. Camera Calibration (K)
fprintf('\n=== Computing Calibration Matrix K (From Normalized Data) ===\n');

if isempty(v_vert_n) || isempty(v_axis_n) || isempty(v_trans_n)
    error('Cannot compute K: missing normalized vanishing points.');
end

v1 = v_vert_n(:);
v2 = v_axis_n(:);
v3 = v_trans_n(:);

fprintf('\n--- DEBUG: Normalized Vanishing Points ---\n');
fprintf('v_vert_n  = [%.4f, %.4f, %.4f]\n', v_vert_n);
fprintf('v_axis_n  = [%.4f, %.4f, %.4f]\n', v_axis_n);
fprintf('v_trans_n = [%.4f, %.4f, %.4f]\n', v_trans_n);

fprintf('\n--- DEBUG: Raw Dot Products (Normalized) ---\n');
fprintf('v1 . v2 = %.4f\n', dot(v1, v2));
fprintf('v1 . v3 = %.4f\n', dot(v1, v3));
fprintf('v2 . v3 = %.4f\n', dot(v2, v3));

build_constraint = @(vi, vj) [vi(1)*vj(1), ...
    vi(2)*vj(2), ...
    vi(1)*vj(3)+vi(3)*vj(1), ...
    vi(2)*vj(3)+vi(3)*vj(2)];

A = zeros(4, 4);
b = zeros(4, 1);

A(1,:) = build_constraint(v1, v2);
b(1) = -v1(3)*v2(3);

A(2,:) = build_constraint(v1, v3);
b(2) = -v1(3)*v3(3);

A(3,:) = build_constraint(v2, v3);
b(3) = -v2(3)*v3(3);

A(4,:) = [1, -1, 0, 0];
b(4) = 0;

fprintf('\n--- DEBUG: Constraint Matrix A (Normalized) ---\n');
fprintf('Condition number of A: %.2e\n', cond(A));

if cond(A) > 1e10
    warning('Matrix A is ill-conditioned (cond = %.2e).', cond(A));
end

x = A \ b;
w1 = x(1); w3 = x(2); w4 = x(3); w5 = x(4); w6 = 1;

omega_n = [w1,  0, w4;
    0, w3, w5;
    w4, w5, w6];

fprintf('\nomega_n (IAC) =\n');
disp(omega_n);

eig_omega = eig(omega_n);
fprintf('Eigenvalues of omega_n: [%.4f, %.4f, %.4f]\n', eig_omega);

if any(eig_omega <= 0)
    warning('omega_n is not positive definite! Relaxing square pixel assumption...');
    A_alt = A(1:3, :);
    b_alt = b(1:3);
    x_alt = pinv(A_alt) * b_alt;
    w1 = x_alt(1); w3 = x_alt(2); w4 = x_alt(3); w5 = x_alt(4);
    omega_n = [w1,  0, w4; 0, w3, w5; w4, w5, w6];
    fprintf('New eigenvalues: [%.4f, %.4f, %.4f]\n', eig(omega_n));
end

try
    R = chol(omega_n);
    K_norm = inv(R);
    K_norm = K_norm / K_norm(3,3);
catch ME
    fprintf('Using eigenvalue decomposition...\n');
    KKt = inv(omega_n);
    KKt = (KKt + KKt') / 2;
    [V, D] = eig(KKt);
    D(D < 0) = abs(D(D < 0));
    K_norm = V * sqrt(D);
    K_norm = K_norm / K_norm(3,3);
end

fprintf('\nK_norm (Normalized K) - Should be Upper Triangular:\n');
disp(K_norm);

fprintf('Denormalizing K (Applying Shift and Scale)...\n');
K = T_n2p * K_norm;
K = K / K(3,3);

fprintf('\nCalibration Matrix K (Pixel Coords):\n');
disp(K);

fx = K(1,1);
fy = K(2,2);
skew = K(1,2);
u0 = K(1,3);
v0 = K(2,3);

fprintf('Intrinsic Parameters:\n');
fprintf('  fx = %.2f px\n', fx);
fprintf('  fy = %.2f px\n', fy);
fprintf('  skew = %.4f\n', skew);
fprintf('  Principal point = (%.2f, %.2f)\n', u0, v0);
fprintf('  Image center = (%.2f, %.2f)\n', cols/2, rows/2);

if fx < 0 || fy < 0
    warning('Negative focal length detected!');
end
if u0 < 0 || u0 > cols || v0 < 0 || v0 > rows
    warning('Principal point outside image bounds!');
end

%% 5. Metric Rectification
fprintf('\n=== Computing Metric Rectification ===\n');

if isempty(l_inf_perp)
    error('Cannot rectify: vanishing line not computed.');
end

% Step 1: Affine rectification - map vanishing line to infinity
l_norm = l_inf_perp(:) / l_inf_perp(3);
H_aff = [1, 0, 0; 0, 1, 0; l_norm(1), l_norm(2), 1];

% Step 2: Metric rectification - use the known perpendicular directions
% Transform the vanishing points through H_aff to get their directions at infinity
v_vert_aff = H_aff * v_vert(:);
v_trans_aff = H_aff * v_trans(:);

% These are now points at infinity (3rd coord ~0). Their direction is (x, y).
dir_vert = v_vert_aff(1:2) / norm(v_vert_aff(1:2));
dir_trans = v_trans_aff(1:2) / norm(v_trans_aff(1:2));

% We want to find S such that:
%   S * dir_vert = [0; 1]  (vertical lines go up)
%   S * dir_trans = [1; 0] (transversal lines go right)
% This is: S * [dir_vert, dir_trans] = [0, 1; 1, 0]
% So: S = [0, 1; 1, 0] * inv([dir_vert, dir_trans])

M = [dir_vert, dir_trans];
S_metric = [0, 1; 1, 0] * inv(M);

% Build the full metric correction homography
H_metric = [S_metric, [0; 0]; 0, 0, 1];

% Combined homography: H_rect = H_metric * H_aff
H_rect = H_metric * H_aff;
fprintf('Applied metric correction: vert->vertical, trans->horizontal.\n');

% Determine a robust output view by sampling points across the image
[xx, yy] = meshgrid(linspace(1, cols, 40), linspace(1, rows, 40));
pts_grid = [xx(:), yy(:), ones(numel(xx), 1)]';

% Use H_aff row to identify points near vanishing line (H_rect may have different structure now)
grid_vals = H_aff(3,:) * pts_grid;

% Pick the dominant side of the vanishing line
side = sign(mean(grid_vals));
if side == 0, side = 1; end

% Filter: points on the correct side and not too close to the line
v_range = abs(max(grid_vals) - min(grid_vals));
margin = 0.15 * v_range;
mask = (sign(grid_vals) == side) & (abs(grid_vals) > margin);

if sum(mask) < 10
    mask = (sign(grid_vals) == side);
end

% Transform grid points to rectified space for limit calculation
pts_trans = H_rect * pts_grid(:, mask);
pts_trans = pts_trans(1:2, :) ./ pts_trans(3, :);

% Find bounding box of visible points in the rotated space
min_x = min(pts_trans(1,:)); max_x = max(pts_trans(1,:));
min_y = min(pts_trans(2,:)); max_y = max(pts_trans(2,:));
w_rect = max_x - min_x; h_rect_val = max_y - min_y;

% Limit extreme aspect ratios which can happen near the vanishing line
if h_rect_val > 5 * w_rect, h_rect_val = 5 * w_rect; end
if w_rect > 5 * h_rect_val, w_rect = 5 * h_rect_val; end

% Final scaling and translation (mapping min_x, min_y to 1, 1)
target_w = 2000;
scale = target_w / w_rect;
T = [scale, 0, -scale*min_x + 1; 0, scale, -scale*min_y + 1; 0, 0, 1];
H_final = T * H_rect;

% Warp and show
out_size = [ceil(scale * h_rect_val), target_w];
fprintf('Rectifying view. Points used for limits: %d. Output size: %d x %d\n', sum(mask), out_size(2), out_size(1));
img_rect = imwarp(img, projective2d(H_final'), 'OutputView', imref2d(out_size));

figure(3); clf;
imshow(img_rect);
title(sprintf('Metric Rectified Image (%dx%d)', out_size(2), out_size(1)), 'FontSize', 14);
hold on;
axis([1, out_size(2), 1, out_size(1)]); % Lock axis to output view

% Helper: check if line is at least partially inside the view
is_in_view = @(p1, p2) all([p1(:); p2(:)] > -500 & [p1(:); p2(:)] < [out_size(2); out_size(1); out_size(2); out_size(1)] + 500);

% Overlay vertical lines in black
if ~isempty(lines_v)
    for i = 1:length(lines_v)
        p1 = H_final*[lines_v(i).p1 1]'; p1 = p1(1:2)/p1(3);
        p2 = H_final*[lines_v(i).p2 1]'; p2 = p2(1:2)/p2(3);
        if is_in_view(p1, p2)
            plot([p1(1) p2(1)], [p1(2) p2(2)], 'k-', 'LineWidth', 2.5);
        end
    end
end

% Overlay transversal lines in white with black outline for visibility
if ~isempty(lines_trans)
    for i = 1:length(lines_trans)
        p1 = H_final*[lines_trans(i).p1 1]'; p1 = p1(1:2)/p1(3);
        p2 = H_final*[lines_trans(i).p2 1]'; p2 = p2(1:2)/p2(3);
        if is_in_view(p1, p2)
            plot([p1(1) p2(1)], [p1(2) p2(2)], 'k-', 'LineWidth', 4); % black outline
            plot([p1(1) p2(1)], [p1(2) p2(2)], 'w-', 'LineWidth', 2); % white core
        end
    end
end

H_R = H_final; % For later sections
fprintf('Rectification complete. Output size: %d x %d\n', out_size(2), out_size(1));

%% 6. 3D Reconstruction
fprintf('\n=== Performing 3D Reconstruction ===\n');

K_inv = inv(K);

if ~isempty(v_vert)
    vert_dir = K_inv * v_vert(:);
    vert_dir = vert_dir / norm(vert_dir);
else
    vert_dir = [0; 1; 0];
end

if ~isempty(v_axis)
    axis_dir = K_inv * v_axis(:);
    axis_dir = axis_dir / norm(axis_dir);
else
    axis_dir = [1; 0; 0];
end

if ~isempty(v_trans)
    trans_dir = K_inv * v_trans(:);
    trans_dir = trans_dir / norm(trans_dir);
else
    trans_dir = [0; 0; 1];
end

fprintf('Direction vectors (camera coords):\n');
fprintf('  vert_dir  = [%.4f, %.4f, %.4f]\n', vert_dir);
fprintf('  axis_dir  = [%.4f, %.4f, %.4f]\n', axis_dir);
fprintf('  trans_dir = [%.4f, %.4f, %.4f]\n', trans_dir);

points_3D = struct('pts', {}, 'arc', {});

for i = 1:length(arcs_A)
    arc_pts_2d = arcs_A{i};
    n_pts = size(arc_pts_2d, 1);
    pts_3d = zeros(n_pts, 3);

    for j = 1:n_pts
        p_img = [arc_pts_2d(j, :), 1]';
        ray = K_inv * p_img;
        ray = ray / norm(ray);
        depth = 10 + i;
        pts_3d(j, :) = (depth * ray)';
    end

    points_3D(end+1).pts = pts_3d;
    points_3D(end).arc = 'A';
end

for i = 1:length(arcs_B)
    arc_pts_2d = arcs_B{i};
    n_pts = size(arc_pts_2d, 1);
    pts_3d = zeros(n_pts, 3);

    for j = 1:n_pts
        p_img = [arc_pts_2d(j, :), 1]';
        ray = K_inv * p_img;
        ray = ray / norm(ray);
        depth = 10 + i;
        pts_3d(j, :) = (depth * ray)';
    end

    points_3D(end+1).pts = pts_3d;
    points_3D(end).arc = 'B';
end

fprintf('Reconstructed %d arcs in 3D\n', length(points_3D));
fprintf('Note: Depths are in arbitrary units. Full reconstruction requires:\n');
fprintf('  - Cylinder fitting to arc families\n');
fprintf('  - Scale constraint from known measurements\n');
fprintf('  - Symmetry constraint between arc families A and B\n');

%% 7. Plotting
figure('Name', '3D Reconstruction'); clf;
hold on;

colors_A = lines(length(arcs_A));
colors_B = lines(length(arcs_B));

for arc_idx = 1:length(points_3D)
    pts = points_3D(arc_idx).pts;

    if strcmp(points_3D(arc_idx).arc, 'A')
        plot3(pts(:,1), pts(:,2), pts(:,3), '.-', 'Color', [0 0.8 0.8], ...
            'LineWidth', 1.5, 'MarkerSize', 10);
    else
        plot3(pts(:,1), pts(:,2), pts(:,3), '.-', 'Color', [0.8 0 0.8], ...
            'LineWidth', 1.5, 'MarkerSize', 10);
    end
end

quiver3(0, 0, 0, axis_dir(1), axis_dir(2), axis_dir(3), 0.3, 'g', 'LineWidth', 2);
quiver3(0, 0, 0, vert_dir(1), vert_dir(2), vert_dir(3), 0.3, 'k', 'LineWidth', 2);
quiver3(0, 0, 0, trans_dir(1), trans_dir(2), trans_dir(3), 0.3, 'r', 'LineWidth', 2);

legend('Arcs A (cyan)', 'Arcs B (magenta)', 'Axis dir', 'Vert dir', 'Trans dir', ...
    'Location', 'best');

grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Reconstruction of Vault Arcs');
view(3);

figure('Name', '3D Views'); clf;
subplot(2,2,1);
for arc_idx = 1:length(points_3D)
    pts = points_3D(arc_idx).pts;
    if strcmp(points_3D(arc_idx).arc, 'A')
        plot3(pts(:,1), pts(:,2), pts(:,3), 'c.-'); hold on;
    else
        plot3(pts(:,1), pts(:,2), pts(:,3), 'm.-'); hold on;
    end
end
view(0, 90); title('Top View (XY)'); axis equal; grid on;

subplot(2,2,2);
for arc_idx = 1:length(points_3D)
    pts = points_3D(arc_idx).pts;
    if strcmp(points_3D(arc_idx).arc, 'A')
        plot3(pts(:,1), pts(:,2), pts(:,3), 'c.-'); hold on;
    else
        plot3(pts(:,1), pts(:,2), pts(:,3), 'm.-'); hold on;
    end
end
view(0, 0); title('Front View (XZ)'); axis equal; grid on;

subplot(2,2,3);
for arc_idx = 1:length(points_3D)
    pts = points_3D(arc_idx).pts;
    if strcmp(points_3D(arc_idx).arc, 'A')
        plot3(pts(:,1), pts(:,2), pts(:,3), 'c.-'); hold on;
    else
        plot3(pts(:,1), pts(:,2), pts(:,3), 'm.-'); hold on;
    end
end
view(90, 0); title('Side View (YZ)'); axis equal; grid on;

subplot(2,2,4);
for arc_idx = 1:length(points_3D)
    pts = points_3D(arc_idx).pts;
    if strcmp(points_3D(arc_idx).arc, 'A')
        plot3(pts(:,1), pts(:,2), pts(:,3), 'c.-'); hold on;
    else
        plot3(pts(:,1), pts(:,2), pts(:,3), 'm.-'); hold on;
    end
end
view(45, 30); title('Isometric View'); axis equal; grid on;

fprintf('\n=== Processing Complete ===\n');
fprintf('Figures generated:\n');
fprintf('  1. All Features Loaded\n');
fprintf('  2. Vanishing Points & Line Extensions\n');
fprintf('  3. Metric Rectified Image\n');
fprintf('  4. 3D Reconstruction of Vault Arcs\n');
fprintf('  5. 3D Views (Top/Front/Side/Isometric)\n');

%% Helper Functions

function lines = select_lines(img, promptStr, colorCode)
if nargin < 3, colorCode = 'y'; end
hFig = figure('Name', promptStr); imshow(img); hold on;
title({promptStr, 'MOUSE: Zoom/Pan with Toolbar tools', ...
    'ACTION: Toggle tools OFF to click/pick points', ...
    '[u]: Undo, [Enter]: Finish category'});

lines = struct('p1', {}, 'p2', {});
plotHandles = {};
tempHandle = [];

while true
    k = waitforbuttonpress;

    if k == 0
        zoomObj = zoom(hFig);
        panObj = pan(hFig);
        if strcmp(zoomObj.Enable, 'on') || strcmp(panObj.Enable, 'on')
            continue;
        end

        currPt = get(gca, 'CurrentPoint');
        x = currPt(1,1);
        y = currPt(1,2);

        if isempty(tempHandle)
            tempHandle = plot(x, y, [colorCode '+'], 'MarkerSize', 12, 'LineWidth', 1.5);
            p1 = [x, y];
            fprintf('Point 1 set at (%.2f, %.2f)\n', x, y);
        else
            p2 = [x, y];
            idx = length(lines) + 1;
            lines(idx).p1 = p1;
            lines(idx).p2 = p2;
            delete(tempHandle);
            tempHandle = [];
            hLine = plot([p1(1), p2(1)], [p1(2), p2(2)], [colorCode '-'], 'LineWidth', 2);
            hP1 = plot(p1(1), p1(2), [colorCode '+'], 'MarkerSize', 8);
            hP2 = plot(p2(1), p2(2), [colorCode '+'], 'MarkerSize', 8);
            plotHandles{end+1} = [hLine, hP1, hP2];
            fprintf('Line added.\n');
        end

    elseif k == 1
        key = get(hFig, 'CurrentCharacter');
        if isempty(key), continue; end

        if key == 13
            if isempty(tempHandle)
                break;
            else
                fprintf('Current line incomplete. Pick second point or [u] to cancel.\n');
            end
        elseif key == 'u'
            if ~isempty(tempHandle)
                delete(tempHandle);
                tempHandle = [];
                fprintf('Point 1 cancelled.\n');
            elseif ~isempty(lines)
                lines(end) = [];
                delete(plotHandles{end});
                plotHandles(end) = [];
                fprintf('Last line undone.\n');
            else
                fprintf('Nothing to undo.\n');
            end
        end
    end
end
close(hFig);
end

function arcs = select_arcs(img, promptStr)
hFig = figure('Name', promptStr); imshow(img); hold on;
title({promptStr, 'MOUSE: Click to pick (ensure Zoom tool is OFF)', ...
    '[Backspace]: Undo point, [u]: Undo Arc, [Enter]: Save Arc/Finish'});

arcs = {};
arcPlots = {};
currentArcPts = [];
currentArcPlot = [];

while true
    title({promptStr, sprintf('Arc %d: %d points picked', length(arcs)+1, size(currentArcPts, 1))});

    k = waitforbuttonpress;
    if k == 0
        zoomObj = zoom(hFig); panObj = pan(hFig);
        if strcmp(zoomObj.Enable, 'on') || strcmp(panObj.Enable, 'on'), continue; end

        currPt = get(gca, 'CurrentPoint');
        x = currPt(1,1); y = currPt(1,2);
        currentArcPts = [currentArcPts; x, y];
        if ~isempty(currentArcPlot), delete(currentArcPlot); end
        currentArcPlot = plot(currentArcPts(:,1), currentArcPts(:,2), 'r.-', 'MarkerSize', 12, 'LineWidth', 1.5);
        fprintf('Point picked at (%.1f, %.1f)\n', x, y);

    elseif k == 1
        key = get(hFig, 'CurrentCharacter');
        if isempty(key), continue; end

        if key == 13
            if isempty(currentArcPts)
                break;
            else
                arcs{end+1} = currentArcPts;
                arcPlots{end+1} = currentArcPlot;
                currentArcPts = [];
                currentArcPlot = [];
                fprintf('Arc saved.\n');
            end
        elseif key == 'u' || key == 8
            if ~isempty(currentArcPts)
                currentArcPts(end,:) = [];
                delete(currentArcPlot);
                if ~isempty(currentArcPts)
                    currentArcPlot = plot(currentArcPts(:,1), currentArcPts(:,2), 'r.-', 'MarkerSize', 12);
                else
                    currentArcPlot = [];
                end
                fprintf('Point undone.\n');
            elseif ~isempty(arcs)
                delete(arcPlots{end});
                arcs(end) = [];
                arcPlots(end) = [];
                fprintf('Whole arc removed.\n');
            end
        end
    end
end
close(hFig);
end

function vp = compute_vanishing_point(lines_hom)
n = length(lines_hom);
if n < 2
    error('Need at least 2 lines to compute a vanishing point.');
end

A = zeros(n, 3);
for k = 1:n
    A(k, :) = lines_hom{k}';
end

[~, ~, V] = svd(A);
vp = V(:, end);

if abs(vp(3)) > 1e-9
    vp = vp / vp(3);
end
vp = vp';
end
%% Verification Script
fprintf('\n=== Verification: Checking Rectified Angles ===\n');

% 1. Check Vertical Lines (Should be 90 degrees wrt X-axis)
angles_v = [];
if ~isempty(lines_v)
    for i = 1:length(lines_v)
        % Map points to rectified space
        p1 = H_final * [lines_v(i).p1 1]'; p1 = p1(1:2)/p1(3);
        p2 = H_final * [lines_v(i).p2 1]'; p2 = p2(1:2)/p2(3);

        % Calculate angle relative to horizontal
        dx = p2(1) - p1(1);
        dy = p2(2) - p1(2);
        angle_deg = atan2d(dy, dx);

        % We expect 90 or -90 (or 270). Take absolute deviation from 90.
        dev = abs(abs(angle_deg) - 90);
        angles_v(end+1) = dev;
    end
    fprintf('Average deviation of Vertical lines from true vertical: %.2f degrees\n', mean(angles_v));
end

% 2. Check Transversal Lines (Should be 0 degrees wrt X-axis)
angles_t = [];
if ~isempty(lines_trans)
    for i = 1:length(lines_trans)
        % Map points to rectified space
        p1 = H_final * [lines_trans(i).p1 1]'; p1 = p1(1:2)/p1(3);
        p2 = H_final * [lines_trans(i).p2 1]'; p2 = p2(1:2)/p2(3);

        % Calculate angle relative to horizontal
        dx = p2(1) - p1(1);
        dy = p2(2) - p1(2);
        angle_deg = atan2d(dy, dx);

        % We expect 0 or 180.
        dev = min(abs(angle_deg), abs(abs(angle_deg)-180));
        angles_t(end+1) = dev;
    end
    fprintf('Average deviation of Transversal lines from true horizontal: %.2f degrees\n', mean(angles_t));
end

if mean([angles_v, angles_t]) < 2.0
    fprintf('>> SUCCESS: Rectification looks geometrically valid.\n');
else
    fprintf('>> WARNING: Rectification might be skewed. Check V_trans and V_vert selection.\n');
end
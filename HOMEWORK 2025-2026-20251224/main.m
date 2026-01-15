% IACV Homework 2025-2026
% 3D Reconstruction of a Cylindric Vault
% Author: [Your Name]
% Date: 2026-01-03

clear; close all; clc;

%% Setup and Parameters
imageFile = 'San Maurizio.jpg';
featureFile = 'features_fixed.mat';

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

if exist(featureFile, 'file')
    fprintf('Loading existing features from %s...\n', featureFile);
    load(featureFile);
else
    fprintf('No feature file found. Starting fresh.\n');
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

% Test all possible vanishing line combinations
fprintf('\n--- Testing Vanishing Line Candidates ---\n');

% Option 1: cross(v_vert, v_trans) - for plane perp to axis (contains vert and trans)
if ~isempty(v_vert) && ~isempty(v_trans)
    l_vt = cross(v_vert, v_trans);
    l_vt = l_vt / norm(l_vt(1:2));
    cv_vt = [l_vt*[1;1;1], l_vt*[cols;1;1], l_vt*[cols;rows;1], l_vt*[1;rows;1]];
    ok_vt = all(cv_vt > 0) || all(cv_vt < 0);
    fprintf('cross(v_vert, v_trans): corners same side = %d\n', ok_vt);
end

% Option 2: cross(v_vert, v_axis) - for plane containing vert and axis
if ~isempty(v_vert) && ~isempty(v_axis)
    l_va = cross(v_vert, v_axis);
    l_va = l_va / norm(l_va(1:2));
    cv_va = [l_va*[1;1;1], l_va*[cols;1;1], l_va*[cols;rows;1], l_va*[1;rows;1]];
    ok_va = all(cv_va > 0) || all(cv_va < 0);
    fprintf('cross(v_vert, v_axis): corners same side = %d\n', ok_va);
end

% Option 3: cross(v_trans, v_axis) - for plane containing trans and axis
if ~isempty(v_trans) && ~isempty(v_axis)
    l_ta = cross(v_trans, v_axis);
    l_ta = l_ta / norm(l_ta(1:2));
    cv_ta = [l_ta*[1;1;1], l_ta*[cols;1;1], l_ta*[cols;rows;1], l_ta*[1;rows;1]];
    ok_ta = all(cv_ta > 0) || all(cv_ta < 0);
    fprintf('cross(v_trans, v_axis): corners same side = %d\n', ok_ta);
end

% Select a vanishing line that has all corners on the same side (for valid rectification)
l_inf_perp = [];
l_inf_name = '';

% Try cross(v_vert, v_trans) first - the plane perpendicular to axis
if ~isempty(v_vert) && ~isempty(v_trans)
    l_test = cross(v_vert, v_trans);
    l_test = l_test / norm(l_test(1:2));
    cv = [l_test*[1;1;1], l_test*[cols;1;1], l_test*[cols;rows;1], l_test*[1;rows;1]];
    if all(cv > 0) || all(cv < 0)
        l_inf_perp = l_test;
        l_inf_name = 'cross(v_vert, v_trans)';
    end
end

% If that doesn't work, try cross(v_trans, v_axis)
if isempty(l_inf_perp) && ~isempty(v_trans) && ~isempty(v_axis)
    l_test = cross(v_trans, v_axis);
    l_test = l_test / norm(l_test(1:2));
    cv = [l_test*[1;1;1], l_test*[cols;1;1], l_test*[cols;rows;1], l_test*[1;rows;1]];
    if all(cv > 0) || all(cv < 0)
        l_inf_perp = l_test;
        l_inf_name = 'cross(v_trans, v_axis)';
        fprintf('NOTE: Using cross(v_trans, v_axis) instead - this is NOT the plane perpendicular to axis!\n');
    end
end

% If that doesn't work either, try cross(v_vert, v_axis)
if isempty(l_inf_perp) && ~isempty(v_vert) && ~isempty(v_axis)
    l_test = cross(v_vert, v_axis);
    l_test = l_test / norm(l_test(1:2));
    cv = [l_test*[1;1;1], l_test*[cols;1;1], l_test*[cols;rows;1], l_test*[1;rows;1]];
    if all(cv > 0) || all(cv < 0)
        l_inf_perp = l_test;
        l_inf_name = 'cross(v_vert, v_axis)';
        fprintf('NOTE: Using cross(v_vert, v_axis) instead - this is NOT the plane perpendicular to axis!\n');
    end
end

if ~isempty(l_inf_perp)
    fprintf('\nSelected vanishing line: %s = [%.4f, %.4f, %.4f]\n', l_inf_name, l_inf_perp);
else
    warning('No valid vanishing line found! All combinations cross through the image.');
end

if ~isempty(v_vert) || ~isempty(v_axis) || ~isempty(v_trans)
    figure(2); clf; imshow(img); title('Vanishing Points & Line Extensions'); hold on;
    all_pts = [1, 1; cols, rows];

    if ~isempty(v_vert) && abs(v_vert(3)) > 1e-9
        vp = v_vert(1:2)/v_vert(3); all_pts = [all_pts; vp];
        plot(vp(1), vp(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
        text(vp(1), vp(2), '  v_{vert}', 'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
        for i=1:length(lines_v)
            mid = (lines_v(i).p1 + lines_v(i).p2) / 2;
            plot([mid(1), vp(1)], [mid(2), vp(2)], 'k:', 'LineWidth', 0.5);
            if ~isempty(l_inf_perp)
                L = cross([lines_v(i).p1 1], [lines_v(i).p2 1]);
                pt_inf = cross(L, l_inf_perp);
                if abs(pt_inf(3)) > 1e-9
                    pt_inf = pt_inf(1:2)/pt_inf(3); all_pts = [all_pts; pt_inf];
                    plot([mid(1), pt_inf(1)], [mid(2), pt_inf(2)], 'k:', 'LineWidth', 0.5);
                end
            end
        end
    end

    if ~isempty(v_axis) && abs(v_axis(3)) > 1e-9
        vp = v_axis(1:2)/v_axis(3); all_pts = [all_pts; vp];
        plot(vp(1), vp(2), 'go', 'MarkerSize', 12, 'LineWidth', 2);
        text(vp(1), vp(2), '  v_{axis}', 'Color', 'g', 'FontSize', 10, 'FontWeight', 'bold');
        for i=1:length(lines_axis)
            mid = (lines_axis(i).p1 + lines_axis(i).p2) / 2;
            plot([mid(1), vp(1)], [mid(2), vp(2)], 'g:', 'LineWidth', 0.5);
            if ~isempty(l_inf_perp)
                L = cross([lines_axis(i).p1 1], [lines_axis(i).p2 1]);
                pt_inf = cross(L, l_inf_perp);
                if abs(pt_inf(3)) > 1e-9
                    pt_inf = pt_inf(1:2)/pt_inf(3); all_pts = [all_pts; pt_inf];
                    plot([mid(1), pt_inf(1)], [mid(2), pt_inf(2)], 'g:', 'LineWidth', 0.5);
                end
            end
        end
        if ~isempty(line_apical)
            mid = (line_apical.p1 + line_apical.p2) / 2;
            plot([mid(1), vp(1)], [mid(2), vp(2)], 'y:', 'LineWidth', 0.8);
            if ~isempty(l_inf_perp)
                L = cross([line_apical.p1 1], [line_apical.p2 1]);
                pt_inf = cross(L, l_inf_perp);
                if abs(pt_inf(3)) > 1e-9
                    pt_inf = pt_inf(1:2)/pt_inf(3); all_pts = [all_pts; pt_inf];
                    plot([mid(1), pt_inf(1)], [mid(2), pt_inf(2)], 'y:', 'LineWidth', 0.8);
                end
            end
        end
    end

    if ~isempty(v_trans) && abs(v_trans(3)) > 1e-9
        vp = v_trans(1:2)/v_trans(3); all_pts = [all_pts; vp];
        plot(vp(1), vp(2), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
        text(vp(1), vp(2), '  v_{trans}', 'Color', [0.8 0 0], 'FontSize', 10, 'FontWeight', 'bold');
        for i=1:length(lines_trans)
            mid = (lines_trans(i).p1 + lines_trans(i).p2) / 2;
            plot([mid(1), vp(1)], [mid(2), vp(2)], 'w:', 'LineWidth', 0.5);
            if ~isempty(l_inf_perp)
                L = cross([lines_trans(i).p1 1], [lines_trans(i).p2 1]);
                pt_inf = cross(L, l_inf_perp);
                if abs(pt_inf(3)) > 1e-9
                    pt_inf = pt_inf(1:2)/pt_inf(3); all_pts = [all_pts; pt_inf];
                    plot([mid(1), pt_inf(1)], [mid(2), pt_inf(2)], 'w:', 'LineWidth', 0.5);
                end
            end
        end
    end

    x_min = min(all_pts(:,1)); x_max = max(all_pts(:,1));
    y_min = min(all_pts(:,2)); y_max = max(all_pts(:,2));
    pad_x = abs(x_max - x_min) * 0.1 + 100;
    pad_y = abs(y_max - y_min) * 0.1 + 100;
    curr_axis = [x_min-pad_x, x_max+pad_x, y_min-pad_y, y_max+pad_y];

    set(gca, 'Visible', 'on');
    axis(curr_axis);
    set(gca, 'Layer', 'top');
    grid on; grid minor;
    set(gca, 'Color', [0.9 0.9 0.9]);

    if ~isempty(l_inf_perp)
        a = l_inf_perp(1); b = l_inf_perp(2); c = l_inf_perp(3);
        xl = [curr_axis(1), curr_axis(2)];
        if abs(b) > 1e-9
            yl = -(a*xl + c)/b;
            plot(xl, yl, 'b-', 'LineWidth', 2.5);
            text(xl(2), yl(2), '  l_{\infty\perp}', 'Color', 'b', 'FontSize', 14, 'FontWeight', 'bold');
        else
            x_vert = -c/a;
            plot([x_vert, x_vert], [curr_axis(3), curr_axis(4)], 'b-', 'LineWidth', 2.5);
            text(x_vert, curr_axis(4), '  l_{\infty\perp}', 'Color', 'b', 'FontSize', 14, 'FontWeight', 'bold');
        end
    end
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
    warning('Cholesky failed on K_norm: %s', ME.message);
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

%% 5. Rectification
fprintf('\n=== Computing Rectification ===\n');

if isempty(l_inf_perp)
    error('Cannot rectify: l_inf_perp not computed.');
end

corners = [1, cols, cols, 1; 1, 1, rows, rows; 1, 1, 1, 1];

% --- Stage 1: Affine Rectification ---
% H_aff maps the vanishing line l_inf_perp to the line at infinity [0; 0; 1]
l = l_inf_perp(:) / l_inf_perp(3);
H_aff = [1, 0, 0; 0, 1, 0; l(1), l(2), 1];

% Affine Safe Warp (for visualization)
c_aff = H_aff * corners; c_aff = c_aff ./ c_aff(3,:);
min_ax = min(c_aff(1,:)); max_ax = max(c_aff(1,:));
min_ay = min(c_aff(2,:)); max_ay = max(c_aff(2,:));
s_aff = 1500 / (max_ax - min_ax);
T_aff = [s_aff, 0, -s_aff*min_ax+1; 0, s_aff, -s_aff*min_ay+1; 0, 0, 1];
H_aff_v = T_aff * H_aff;
img_aff = imwarp(img, projective2d(H_aff_v'), 'OutputView', imref2d([ceil(s_aff*(max_ay-min_ay)), 1500]));

figure(3); clf; imshow(img_aff); title('Step 1: Affine Rectified Image');

% --- Stage 2: Euclidean Rectification (H_R) ---
% For a plane perpendicular to the cylinder axis, containing vertical and transversal directions
% The vanishing line of this plane is l_inf_perp = cross(v_vert, v_trans)

fprintf('\n=== Euclidean Rectification ===\n');

% Check if vanishing line passes through the image
l = l_inf_perp(:);
corner_vals = [l'*[1;1;1], l'*[cols;1;1], l'*[cols;rows;1], l'*[1;rows;1]];
fprintf('Corner values relative to vanishing line: [%.1f, %.1f, %.1f, %.1f]\n', corner_vals);

% Determine which side of the vanishing line has the most image content
positive_side = sum(corner_vals > 0);
negative_side = sum(corner_vals < 0);
fprintf('Corners on positive side: %d, negative side: %d\n', positive_side, negative_side);

if positive_side == 4 || negative_side == 4
    fprintf('All corners on same side - full image can be rectified.\n');
    % Use the full image
    rect_corners = corners;
    use_full_image = true;
else
    warning('Vanishing line passes through image! This is a DEGENERATE case.');
    fprintf('The vanishing line l_inf_perp crosses through the image.\n');
    fprintf('This means the plane being rectified is viewed almost edge-on.\n');
    fprintf('Full perspective rectification is not possible.\n\n');

    % For visualization, we'll just show the affine rectified image (Figure 3)
    % and skip the problematic full rectification
    use_full_image = false;
end

if use_full_image
    % Standard affine rectification: map vanishing line to infinity
    l_norm = l_inf_perp(:) / l_inf_perp(3);
    H_R = [1, 0, 0; 0, 1, 0; l_norm(1), l_norm(2), 1];

    % Transform corners
    pts_out = H_R * corners;
    pts_out = pts_out ./ pts_out(3,:);

    % Auto-rotate to make vertical lines vertical
    if ~isempty(lines_v)
        p1 = [lines_v(1).p1, 1]';
        p2 = [lines_v(1).p2, 1]';
        q1 = H_R * p1; q1 = q1(1:2)/q1(3);
        q2 = H_R * p2; q2 = q2(1:2)/q2(3);

        if q2(2) > q1(2)
            dx = q1(1) - q2(1); dy = q1(2) - q2(2);
        else
            dx = q2(1) - q1(1); dy = q2(2) - q1(2);
        end

        theta = atan2(dy, dx);
        rot_angle = -pi/2 - theta;
        while rot_angle > pi, rot_angle = rot_angle - 2*pi; end
        while rot_angle < -pi, rot_angle = rot_angle + 2*pi; end

        c_r = cos(rot_angle); s_r = sin(rot_angle);
        H_rot = [c_r, -s_r, 0; s_r, c_r, 0; 0, 0, 1];
        H_R = H_rot * H_R;
        pts_out = H_R * corners;
        pts_out = pts_out ./ pts_out(3,:);
        fprintf('Auto-rotation: %.1f deg\n', rad2deg(rot_angle));
    end

    % Scale and translate to fit
    min_x = min(pts_out(1,:)); max_x = max(pts_out(1,:));
    min_y = min(pts_out(2,:)); max_y = max(pts_out(2,:));
    w_raw = max_x - min_x;
    h_raw = max_y - min_y;

    target_size = 2000;
    scale_factor = target_size / max(w_raw, h_raw);
    T_final = [scale_factor, 0, -scale_factor*min_x;
        0, scale_factor, -scale_factor*min_y;
        0, 0, 1];
    H_R = T_final * H_R;

    out_width = ceil(scale_factor * w_raw) + 10;
    out_height = ceil(scale_factor * h_raw) + 10;

    fprintf('Output size: %d x %d\n', out_width, out_height);

    outputView = imref2d([out_height, out_width]);
    img_rect = imwarp(img, projective2d(H_R'), 'OutputView', outputView);

    figure(4); clf;
    imshow(img_rect);
    title('Euclidean Rectified Image');
    hold on;
else
    % Degenerate case: just copy the affine result
    fprintf('\nUsing affine rectification result (Figure 3) instead.\n');
    fprintf('Note: Full Euclidean rectification not possible for this view.\n');

    % Display the same as Figure 3 for Figure 4
    H_R = H_aff_v;  % Use the affine homography with scaling
    img_rect = img_aff;

    figure(4); clf;
    imshow(img_rect);
    title('Affine Rectified (Euclidean not possible - vanishing line in image)');
    hold on;
end

fprintf('H_R = \n');
disp(H_R);

if ~isempty(lines_v)
    for i = 1:length(lines_v)
        p1_h = [lines_v(i).p1, 1]';
        p2_h = [lines_v(i).p2, 1]';
        p1_r = H_R * p1_h; p1_r = p1_r(1:2)'/p1_r(3);
        p2_r = H_R * p2_h; p2_r = p2_r(1:2)'/p2_r(3);
        plot([p1_r(1) p2_r(1)], [p1_r(2) p2_r(2)], 'k-', 'LineWidth', 2);
    end
end
if ~isempty(lines_trans)
    for i = 1:length(lines_trans)
        p1_h = [lines_trans(i).p1, 1]';
        p2_h = [lines_trans(i).p2, 1]';
        p1_r = H_R * p1_h; p1_r = p1_r(1:2)'/p1_r(3);
        p2_r = H_R * p2_h; p2_r = p2_r(1:2)'/p2_r(3);
        plot([p1_r(1) p2_r(1)], [p1_r(2) p2_r(2)], 'w-', 'LineWidth', 2);
    end
end

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
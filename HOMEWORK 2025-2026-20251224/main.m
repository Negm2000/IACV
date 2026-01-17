% IACV Homework 2025-2026
% 3D Reconstruction of a Cylindric Vault
% Refactored Main Script
%
% Dependencies: feature_picker.m, compute_vanishing_geometry.m,
%               compute_calibration_matrix.m, compute_metric_rectification.m,
%               reconstruct_3d.m

clear; close all; clc;

%% Configuration
imageFile = 'San Maurizio.jpg';
featureFile = 'features_fixed.mat';

%% 1. Load Image and Features
fprintf('=== Loading Image and Features ===\n');
img = imread(imageFile);
[rows, cols, ~] = size(img);

if ~exist(featureFile, 'file')
    error('Run feature_picker.m first to select features.');
end
load(featureFile);

fprintf('Loaded: %d vertical, %d axis, %d trans lines, %d+%d arcs\n', ...
    length(lines_v), length(lines_axis), length(lines_trans), length(arcs_A), length(arcs_B));

%% 2. Compute Vanishing Geometry
fprintf('\n=== Computing Vanishing Points ===\n');
[v_vert, v_axis, v_trans, v_vert_n, v_axis_n, v_trans_n, l_inf_perp] = ...
    compute_vanishing_geometry(img, lines_v, lines_axis, line_apical, lines_trans);

fprintf('v_vert = [%.2f, %.2f, %.4f]\n', v_vert);
fprintf('v_axis = [%.2f, %.2f, %.4f]\n', v_axis);
fprintf('v_trans = [%.2f, %.2f, %.4f]\n', v_trans);
fprintf('l_inf_perp = [%.4f, %.4f, %.4f]\n', l_inf_perp);

%% 3. Compute Calibration Matrix K
fprintf('\n=== Computing Calibration Matrix ===\n');
K = compute_calibration_matrix(v_vert_n, v_axis_n, v_trans_n, img);
fprintf('K =\n'); disp(K);
fprintf('fx=%.1f, fy=%.1f, (u0,v0)=(%.1f,%.1f)\n', K(1,1), K(2,2), K(1,3), K(2,3));

%% 4. Metric Rectification
fprintf('\n=== Computing Metric Rectification ===\n');
[H_R, img_rect, out_size] = compute_metric_rectification(img, l_inf_perp, v_vert, v_trans);
fprintf('Output size: %dx%d\n', out_size(2), out_size(1));

%% 5. 3D Reconstruction
fprintf('\n=== 3D Reconstruction ===\n');
[points_3D, axis_dir, vert_dir, trans_dir, P_axis] = ...
    reconstruct_3d(K, v_vert, v_axis, v_trans, arcs_A, arcs_B, line_apical);
fprintf('Reconstructed %d arcs, scale calibrated to d=1\n', length(points_3D));

%% 6. Find Nodal Points
fprintf('\n=== Finding Nodal Points ===\n');
nodal_points = find_nodal_points(arcs_A, arcs_B, H_R);
fprintf('Found %d nodal points\n', size(nodal_points, 1));

%% 7. Visualization
plot_all_results(img, img_rect, out_size, H_R, ...
    lines_v, lines_axis, line_apical, lines_trans, arcs_A, arcs_B, ...
    v_vert, v_axis, v_trans, l_inf_perp, ...
    points_3D, axis_dir, vert_dir, trans_dir, P_axis, nodal_points);

%% 8. Verification
verify_rectification(lines_v, lines_trans, H_R);

fprintf('\n=== Complete ===\n');

%% ========== HELPER FUNCTIONS ==========

function nodal_points = find_nodal_points(arcs_A, arcs_B, H_R)
nodal_points = [];
if isempty(arcs_A) || isempty(arcs_B), return; end

for i = 1:length(arcs_A)
    pts_A = arcs_A{i};
    for j = 1:length(arcs_B)
        pts_B = arcs_B{j};
        for k = 1:size(pts_A, 1)
            p_A = H_R * [pts_A(k, :), 1]';
            p_A_rect = p_A(1:2) / p_A(3);
            for m = 1:size(pts_B, 1)
                p_B = H_R * [pts_B(m, :), 1]';
                p_B_rect = p_B(1:2) / p_B(3);
                if norm(p_A_rect - p_B_rect) < 50
                    mid_img = (pts_A(k, :) + pts_B(m, :)) / 2;
                    nodal_points = [nodal_points; mid_img, i, j];
                end
            end
        end
    end
end
if ~isempty(nodal_points)
    [~, idx] = unique(nodal_points(:, 3:4), 'rows', 'first');
    nodal_points = nodal_points(idx, :);
end
end

function plot_all_results(img, img_rect, out_size, H_R, ...
    lines_v, lines_axis, line_apical, lines_trans, arcs_A, arcs_B, ...
    v_vert, v_axis, v_trans, l_inf_perp, ...
    points_3D, axis_dir, vert_dir, trans_dir, P_axis, nodal_points)

[rows, cols, ~] = size(img);

% Figure 1: Features
figure(1); clf; imshow(img); title('Features'); hold on;
plot_features(lines_v, lines_axis, line_apical, lines_trans, arcs_A, arcs_B);

% Figure 2: Vanishing Points
figure(2); clf; imshow(img); title('Vanishing Geometry'); hold on;
plot_features(lines_v, lines_axis, line_apical, lines_trans, arcs_A, arcs_B);
if ~isempty(l_inf_perp), draw_vanishing_line(l_inf_perp, rows, cols); end

% Figure 3: Rectified Image
figure(3); clf; imshow(img_rect);
title(sprintf('Rectified (%dx%d)', out_size(2), out_size(1))); hold on;

% Figure 4: Nodal Points
figure(4); clf; imshow(img_rect); title('Nodal Points'); hold on;
for n = 1:size(nodal_points, 1)
    np = nodal_points(n, :);
    p = H_R * [np(1:2), 1]'; p = p(1:2)/p(3);
    plot(p(1), p(2), 'yo', 'MarkerSize', 12, 'LineWidth', 2);
end

% Figure 5: 3D Reconstruction
figure(5); clf; hold on;
for i = 1:length(points_3D)
    pts = points_3D(i).pts;
    if strcmp(points_3D(i).arc, 'A')
        plot3(pts(:,1), pts(:,2), pts(:,3), 'c.-', 'LineWidth', 1.5);
    else
        plot3(pts(:,1), pts(:,2), pts(:,3), 'm.-', 'LineWidth', 1.5);
    end
end
if ~isempty(P_axis)
    t = linspace(-5, 5, 50)';
    axis_line = P_axis' + t * axis_dir';
    plot3(axis_line(:,1), axis_line(:,2), axis_line(:,3), 'g-', 'LineWidth', 3);
end
quiver3(0,0,0, axis_dir(1), axis_dir(2), axis_dir(3), 2, 'g', 'LineWidth', 2);
quiver3(0,0,0, vert_dir(1), vert_dir(2), vert_dir(3), 2, 'k', 'LineWidth', 2);
quiver3(0,0,0, trans_dir(1), trans_dir(2), trans_dir(3), 2, 'r', 'LineWidth', 2);
grid on; axis equal; xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Reconstruction'); view(3);
legend('Arcs A', 'Arcs B', 'Cylinder Axis', 'axis dir', 'vert dir', 'trans dir');

% Figure 6: Multi-view
figure(6); clf;
views = {[0 90], [0 0], [90 0], [45 30]};
titles = {'Top', 'Front', 'Side', 'Iso'};
for v = 1:4
    subplot(2,2,v); hold on;
    for i = 1:length(points_3D)
        pts = points_3D(i).pts;
        if strcmp(points_3D(i).arc, 'A')
            plot3(pts(:,1), pts(:,2), pts(:,3), 'c.-');
        else
            plot3(pts(:,1), pts(:,2), pts(:,3), 'm.-');
        end
    end
    view(views{v}); title(titles{v}); axis equal; grid on;
end
sgtitle('3D Views');

% Figure 7: Single Arc
if ~isempty(points_3D)
    figure(7); clf;
    pts = points_3D(1).pts;
    for v = 1:4
        subplot(2,2,v);
        plot3(pts(:,1), pts(:,2), pts(:,3), 'c.-', 'LineWidth', 2);
        view(views{v}); title(titles{v}); axis equal; grid on;
        xlabel('X'); ylabel('Y'); zlabel('Z');
    end
    sgtitle('Single Arc (A1)');
    fprintf('Arc A1 3D coords (%d pts):\n', size(pts,1));
    for p = 1:min(12, size(pts,1))
        fprintf('  [%.3f, %.3f, %.3f]\n', pts(p,:));
    end
end
end

function plot_features(lines_v, lines_axis, line_apical, lines_trans, arcs_A, arcs_B)
for i=1:length(lines_v), plot([lines_v(i).p1(1) lines_v(i).p2(1)], [lines_v(i).p1(2) lines_v(i).p2(2)], 'k-', 'LineWidth', 2); end
for i=1:length(lines_axis), plot([lines_axis(i).p1(1) lines_axis(i).p2(1)], [lines_axis(i).p1(2) lines_axis(i).p2(2)], 'g-', 'LineWidth', 2); end
if ~isempty(line_apical), plot([line_apical.p1(1) line_apical.p2(1)], [line_apical.p1(2) line_apical.p2(2)], 'y-', 'LineWidth', 3); end
for i=1:length(lines_trans), plot([lines_trans(i).p1(1) lines_trans(i).p2(1)], [lines_trans(i).p1(2) lines_trans(i).p2(2)], 'w-', 'LineWidth', 2); end
for i=1:length(arcs_A), pts=arcs_A{i}; plot(pts(:,1), pts(:,2), 'c.-', 'MarkerSize', 8); end
for i=1:length(arcs_B), pts=arcs_B{i}; plot(pts(:,1), pts(:,2), 'm.-', 'MarkerSize', 8); end
end

function draw_vanishing_line(l, rows, cols)
a = l(1); b = l(2); c = l(3);
if abs(b) > 1e-6
    x = [-cols*5, cols*6];
    y = (-c - a*x) / b;
    plot(x, y, 'b-', 'LineWidth', 2.5);
else
    x_v = -c/a;
    plot([x_v x_v], [-rows*5, rows*5], 'b-', 'LineWidth', 2.5);
end
end

function verify_rectification(lines_v, lines_trans, H_R)
fprintf('\n=== Verification ===\n');
angles_v = [];
for i = 1:length(lines_v)
    p1 = H_R * [lines_v(i).p1 1]'; p1 = p1(1:2)/p1(3);
    p2 = H_R * [lines_v(i).p2 1]'; p2 = p2(1:2)/p2(3);
    angles_v(end+1) = abs(abs(atan2d(p2(2)-p1(2), p2(1)-p1(1))) - 90);
end
angles_t = [];
for i = 1:length(lines_trans)
    p1 = H_R * [lines_trans(i).p1 1]'; p1 = p1(1:2)/p1(3);
    p2 = H_R * [lines_trans(i).p2 1]'; p2 = p2(1:2)/p2(3);
    angle = atan2d(p2(2)-p1(2), p2(1)-p1(1));
    angles_t(end+1) = min(abs(angle), abs(abs(angle)-180));
end
fprintf('Vertical deviation: %.2f deg, Horizontal deviation: %.2f deg\n', mean(angles_v), mean(angles_t));
if mean([angles_v, angles_t]) < 2
    fprintf('SUCCESS: Rectification valid\n');
else
    fprintf('WARNING: Check feature selection\n');
end
end
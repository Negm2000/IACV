% IACV Homework 2025-2026 - 3D Reconstruction of a Cylindric Vault

clear; close all; clc;

%% Configuration
imageFile = 'San Maurizio.jpg';
featureFile = 'features_fixed.mat';

%% 1. Load Image and Features
fprintf('=== Loading ===\n');
img = imread(imageFile);
[rows, cols, ~] = size(img);
if ~exist(featureFile, 'file'), error('Run feature_picker.m first'); end
load(featureFile);
fprintf('Loaded: %d vert, %d axis, %d trans, %d+%d arcs\n', ...
    length(lines_v), length(lines_axis), length(lines_trans), length(arcs_A), length(arcs_B));

%% 2. Vanishing Geometry
fprintf('\n=== Vanishing Points ===\n');
[v_vert, v_axis, v_trans, v_vert_n, v_axis_n, v_trans_n, l_inf_perp] = ...
    compute_vanishing_geometry(img, lines_v, lines_axis, line_apical, lines_trans);

%% 3. Calibration Matrix K
fprintf('\n=== Calibration ===\n');
K = compute_calibration_matrix(v_vert_n, v_axis_n, v_trans_n, img);
fprintf('K: fx=%.0f, fy=%.0f, (u0,v0)=(%.0f,%.0f)\n', K(1,1), K(2,2), K(1,3), K(2,3));

%% 4. Metric Rectification
fprintf('\n=== Rectification ===\n');
[H_R, img_rect, out_size] = compute_metric_rectification(img, l_inf_perp, v_vert, v_trans);
fprintf('Output: %dx%d\n', out_size(2), out_size(1));

%% 5. 3D Reconstruction
fprintf('\n=== 3D Reconstruction ===\n');
[points_3D, axis_dir, vert_dir, trans_dir, P_axis] = ...
    reconstruct_3d(K, v_vert, v_axis, v_trans, arcs_A, arcs_B, line_apical);
fprintf('%d arcs reconstructed\n', length(points_3D));

%% 6. Nodal Points (Geometric - using vault symmetry)
fprintf('\n=== Nodal Points ===\n');
nodal_points = find_nodal_points_geometric(arcs_A, arcs_B, H_R, v_axis, v_vert);

%% 7. Plots

% Fig 1: Features
figure(1); imshow(img); title('Features'); hold on;
for i=1:length(lines_v), L=lines_v(i); plot([L.p1(1) L.p2(1)], [L.p1(2) L.p2(2)], 'k-', 'LineWidth', 2); end
for i=1:length(lines_axis), L=lines_axis(i); plot([L.p1(1) L.p2(1)], [L.p1(2) L.p2(2)], 'g-', 'LineWidth', 2); end
if ~isempty(line_apical), plot([line_apical.p1(1) line_apical.p2(1)], [line_apical.p1(2) line_apical.p2(2)], 'y-', 'LineWidth', 3); end
for i=1:length(lines_trans), L=lines_trans(i); plot([L.p1(1) L.p2(1)], [L.p1(2) L.p2(2)], 'w-', 'LineWidth', 2); end
for i=1:length(arcs_A), pts=arcs_A{i}; plot(pts(:,1), pts(:,2), 'c.-'); end
for i=1:length(arcs_B), pts=arcs_B{i}; plot(pts(:,1), pts(:,2), 'm.-'); end

% Fig 2: Vanishing Points and Line
figure(2); clf; imshow(img); title('Vanishing Points & Line'); hold on;
% Draw vanishing line
a=l_inf_perp(1); b=l_inf_perp(2); c=l_inf_perp(3);
if abs(b)>1e-6, x=[-cols*2, cols*3]; plot(x, (-c-a*x)/b, 'b-', 'LineWidth', 2); end
% Draw lines converging to vanishing points
if ~isempty(v_vert) && abs(v_vert(3)) > 1e-6
    vp = v_vert(1:2)/v_vert(3);
    for i=1:length(lines_v)
        plot([lines_v(i).p1(1), vp(1)], [lines_v(i).p1(2), vp(2)], 'k:', 'LineWidth', 0.5);
    end
    plot(vp(1), vp(2), 'ko', 'MarkerSize', 15, 'LineWidth', 3);
end
if ~isempty(v_axis) && abs(v_axis(3)) > 1e-6
    vp = v_axis(1:2)/v_axis(3);
    for i=1:length(lines_axis)
        plot([lines_axis(i).p1(1), vp(1)], [lines_axis(i).p1(2), vp(2)], 'g:', 'LineWidth', 0.5);
    end
    plot(vp(1), vp(2), 'go', 'MarkerSize', 15, 'LineWidth', 3);
end
if ~isempty(v_trans) && abs(v_trans(3)) > 1e-6
    vp = v_trans(1:2)/v_trans(3);
    for i=1:length(lines_trans)
        plot([lines_trans(i).p1(1), vp(1)], [lines_trans(i).p1(2), vp(2)], 'r:', 'LineWidth', 0.5);
    end
    plot(vp(1), vp(2), 'ro', 'MarkerSize', 15, 'LineWidth', 3);
end
% Extend axes to show vanishing points
axis auto; ax = axis; axis([min(ax(1),-500) max(ax(2),cols+500) min(ax(3),-500) max(ax(4),rows+500)]);
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.3, 'GridColor', [0.5 0.5 0.5]);
set(gca, 'Layer', 'top');  % Grid on top of image

% Fig 3: Rectified
figure(3); clf; imshow(img_rect); title(sprintf('Rectified %dx%d', out_size(2), out_size(1)));

% Fig 4: Nodal points
figure(4); clf; imshow(img_rect); title('Nodal Points'); hold on;
for n=1:size(nodal_points,1)
    p = H_R * [nodal_points(n,1:2), 1]'; p = p(1:2)/p(3);
    plot(p(1), p(2), 'yo', 'MarkerSize', 12, 'LineWidth', 2);
    text(p(1)+10, p(2), sprintf('N%d%d', nodal_points(n,3), nodal_points(n,4)), 'Color', 'y');
end

% Fig 5: 3D
fig5 = figure(5); clf; hold on;
set(fig5, 'Visible', 'on', 'WindowState', 'normal');
movegui(fig5, 'center');
for i=1:length(points_3D)
    pts = points_3D(i).pts;
    if strcmp(points_3D(i).arc, 'A'), c='c'; else, c='m'; end
    plot3(pts(:,1), pts(:,2), pts(:,3), [c '.-'], 'LineWidth', 1.5);
end
if ~isempty(P_axis)
    t = linspace(-5,5,50)';
    ax = P_axis' + t*axis_dir';
    plot3(ax(:,1), ax(:,2), ax(:,3), 'g-', 'LineWidth', 3);
end
quiver3(0,0,0, axis_dir(1), axis_dir(2), axis_dir(3), 2, 'g', 'LineWidth', 2);
grid on; axis equal; xlabel('X'); ylabel('Y'); zlabel('Z'); title('3D Reconstruction'); view(3);
drawnow;

% Fig 6: Views
fig6 = figure(6); clf;
set(fig6, 'Visible', 'on', 'WindowState', 'normal');
movegui(fig6, 'northeast');
vs = {[0 90], [0 0], [90 0], [45 30]};
for v=1:4
    subplot(2,2,v); hold on;
    for i=1:length(points_3D)
        pts = points_3D(i).pts;
        if points_3D(i).arc=='A', c='c'; else, c='m'; end
        plot3(pts(:,1), pts(:,2), pts(:,3), [c '.-']);
    end
    view(vs{v}); axis equal; grid on;
end
sgtitle('3D Views');

% Fig 7: Single arc
if ~isempty(points_3D)
    fig7 = figure(7); clf;
    set(fig7, 'Visible', 'on', 'WindowState', 'normal');
    movegui(fig7, 'southeast');
    pts = points_3D(1).pts;
    for v=1:4
        subplot(2,2,v);
        plot3(pts(:,1), pts(:,2), pts(:,3), 'c.-', 'LineWidth', 2);
        view(vs{v}); axis equal; grid on;
    end
    sgtitle('Arc A1');
    fprintf('\nArc A1 coords:\n');
    for p=1:min(12,size(pts,1)), fprintf('[%.2f,%.2f,%.2f]\n', pts(p,:)); end
end

%% 8. Verification
angles_v = []; angles_t = [];
for i=1:length(lines_v)
    p1 = H_R*[lines_v(i).p1 1]'; p2 = H_R*[lines_v(i).p2 1]';
    p1=p1(1:2)/p1(3); p2=p2(1:2)/p2(3);
    angles_v(end+1) = abs(abs(atan2d(p2(2)-p1(2), p2(1)-p1(1))) - 90);
end
for i=1:length(lines_trans)
    p1 = H_R*[lines_trans(i).p1 1]'; p2 = H_R*[lines_trans(i).p2 1]';
    p1=p1(1:2)/p1(3); p2=p2(1:2)/p2(3);
    a = atan2d(p2(2)-p1(2), p2(1)-p1(1));
    angles_t(end+1) = min(abs(a), abs(abs(a)-180));
end
fprintf('\n=== Verification ===\n');
fprintf('Vert: %.2f deg, Horiz: %.2f deg\n', mean(angles_v), mean(angles_t));
if mean([angles_v angles_t]) < 2, fprintf('SUCCESS\n'); else, fprintf('WARNING\n'); end

fprintf('\n=== Done ===\n');
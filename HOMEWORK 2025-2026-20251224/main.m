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

%% 6. Nodal Points
nodal_points = find_nodal_points(arcs_A, arcs_B, H_R);
fprintf('%d nodal points\n', size(nodal_points, 1));

%% 7. Plots

% Fig 1: Features
figure(1); imshow(img); title('Features'); hold on;
for i=1:length(lines_v), L=lines_v(i); plot([L.p1(1) L.p2(1)], [L.p1(2) L.p2(2)], 'k-', 'LineWidth', 2); end
for i=1:length(lines_axis), L=lines_axis(i); plot([L.p1(1) L.p2(1)], [L.p1(2) L.p2(2)], 'g-', 'LineWidth', 2); end
if ~isempty(line_apical), plot([line_apical.p1(1) line_apical.p2(1)], [line_apical.p1(2) line_apical.p2(2)], 'y-', 'LineWidth', 3); end
for i=1:length(lines_trans), L=lines_trans(i); plot([L.p1(1) L.p2(1)], [L.p1(2) L.p2(2)], 'w-', 'LineWidth', 2); end
for i=1:length(arcs_A), pts=arcs_A{i}; plot(pts(:,1), pts(:,2), 'c.-'); end
for i=1:length(arcs_B), pts=arcs_B{i}; plot(pts(:,1), pts(:,2), 'm.-'); end

% Fig 2: Vanishing line
figure(2); imshow(img); title('Vanishing Line'); hold on;
a=l_inf_perp(1); b=l_inf_perp(2); c=l_inf_perp(3);
if abs(b)>1e-6, x=[-cols*5, cols*6]; plot(x, (-c-a*x)/b, 'b-', 'LineWidth', 2); end

% Fig 3: Rectified
figure(3); imshow(img_rect); title(sprintf('Rectified %dx%d', out_size(2), out_size(1)));

% Fig 4: Nodal points
figure(4); imshow(img_rect); title('Nodal Points'); hold on;
for n=1:size(nodal_points,1)
    p = H_R * [nodal_points(n,1:2), 1]'; p = p(1:2)/p(3);
    plot(p(1), p(2), 'yo', 'MarkerSize', 12, 'LineWidth', 2);
end

% Fig 5: 3D
figure(5); hold on;
for i=1:length(points_3D)
    pts = points_3D(i).pts;
    if points_3D(i).arc=='A', c='c'; else, c='m'; end
    plot3(pts(:,1), pts(:,2), pts(:,3), [c '.-'], 'LineWidth', 1.5);
end
if ~isempty(P_axis)
    t = linspace(-5,5,50)';
    ax = P_axis' + t*axis_dir';
    plot3(ax(:,1), ax(:,2), ax(:,3), 'g-', 'LineWidth', 3);
end
quiver3(0,0,0, axis_dir(1), axis_dir(2), axis_dir(3), 2, 'g', 'LineWidth', 2);
grid on; axis equal; xlabel('X'); ylabel('Y'); zlabel('Z'); title('3D'); view(3);

% Fig 6: Views
figure(6);
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
    figure(7); pts = points_3D(1).pts;
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

%% Helper
function np = find_nodal_points(arcs_A, arcs_B, H_R)
np = [];
for i=1:length(arcs_A), for j=1:length(arcs_B)
        pA=arcs_A{i}; pB=arcs_B{j};
        for k=1:size(pA,1), for m=1:size(pB,1)
                a = H_R*[pA(k,:) 1]'; a=a(1:2)/a(3);
                b = H_R*[pB(m,:) 1]'; b=b(1:2)/b(3);
                if norm(a-b)<50, np=[np; (pA(k,:)+pB(m,:))/2, i, j]; end
        end, end
end, end
if ~isempty(np), [~,idx]=unique(np(:,3:4),'rows','first'); np=np(idx,:); end
end
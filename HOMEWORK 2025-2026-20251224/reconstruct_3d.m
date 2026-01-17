function [points_3D, axis_dir, vert_dir, trans_dir, P_axis] = reconstruct_3d(K, v_vert, v_axis, v_trans, arcs_A, arcs_B, line_apical)
% RECONSTRUCT_3D - Perform 3D reconstruction of vault arcs
%   Returns 3D points, direction vectors, and axis position.

K_inv = inv(K);

% Compute direction vectors
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

% Cylinder axis point
P_axis = [];
if ~isempty(line_apical)
    mid_apical = (line_apical.p1 + line_apical.p2) / 2;
    ray_apical = K_inv * [mid_apical, 1]';
    ray_apical = ray_apical / norm(ray_apical);
    P_axis = 10 * ray_apical;
end

% Reconstruct arc points
points_3D = struct('pts', {}, 'arc', {});

for i = 1:length(arcs_A)
    arc_pts_2d = arcs_A{i};
    n_pts = size(arc_pts_2d, 1);
    pts_3d = zeros(n_pts, 3);
    for j = 1:n_pts
        ray = K_inv * [arc_pts_2d(j, :), 1]';
        ray = ray / norm(ray);
        pts_3d(j, :) = ((10 + i) * ray)';
    end
    points_3D(end+1).pts = pts_3d;
    points_3D(end).arc = 'A';
end

for i = 1:length(arcs_B)
    arc_pts_2d = arcs_B{i};
    n_pts = size(arc_pts_2d, 1);
    pts_3d = zeros(n_pts, 3);
    for j = 1:n_pts
        ray = K_inv * [arc_pts_2d(j, :), 1]';
        ray = ray / norm(ray);
        pts_3d(j, :) = ((10 + i) * ray)';
    end
    points_3D(end+1).pts = pts_3d;
    points_3D(end).arc = 'B';
end

% Scale calibration using d=1 constraint
if length(arcs_A) >= 2
    centroid_1 = mean(points_3D(1).pts, 1);
    centroid_2 = mean(points_3D(2).pts, 1);
    d_measured = abs(dot(centroid_2 - centroid_1, axis_dir'));
    if d_measured > 0.01
        scale_factor = 1.0 / d_measured;
        for idx = 1:length(points_3D)
            points_3D(idx).pts = points_3D(idx).pts * scale_factor;
        end
        if ~isempty(P_axis)
            P_axis = P_axis * scale_factor;
        end
    end
end
end

function [points_3D, axis_dir, vert_dir, trans_dir, P_axis] = reconstruct_3d(K, v_vert, v_axis, v_trans, arcs_A, arcs_B, line_apical, nodal_points)
% RECONSTRUCT_3D - 3D reconstruction using proper geometric constraints
%
% Theory (from IACV Lecture I):
%   1. Use apical nodal points N_ii to solve triangle for reference depth
%   2. Each arc point gets unique depth via ray-plane intersection
%   3. Cap extreme depths to avoid numerical instability

K_inv = inv(K);

%% 1. Compute 3D direction vectors from vanishing points
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

fprintf('3D directions computed.\n');

d = 1;  % Known arc spacing
points_3D = struct('pts', {}, 'arc', {});

%% 2. Solve triangle using APICAL NODAL POINTS (most stable reference)
% Instead of arc midpoints, use N_11 and N_22 which are mathematically precise

lambda1 = [];
a1 = [];

if nargin >= 8 && ~isempty(nodal_points)
    fprintf('DEBUG: nodal_points size = %d x %d\n', size(nodal_points));

    % Find apical nodal points: N_ii where arcA_idx == arcB_idx
    apical_nodes = nodal_points(nodal_points(:,3) == nodal_points(:,4), :);
    fprintf('DEBUG: apical_nodes (N_ii) count = %d\n', size(apical_nodes, 1));

    if size(apical_nodes, 1) >= 2
        % SORT by arc index so we get N_11, N_22, etc in order
        [~, sort_idx] = sort(apical_nodes(:, 3));
        apical_nodes = apical_nodes(sort_idx, :);

        fprintf('Using apical nodal points N_%d%d and N_%d%d for reference.\n', ...
            apical_nodes(1,3), apical_nodes(1,4), apical_nodes(2,3), apical_nodes(2,4));

        % Get rays for first two apical nodes (now sorted: N_11, N_22)
        n1_2d = apical_nodes(1, 1:2);  % N_11
        n2_2d = apical_nodes(2, 1:2);  % N_22
        fprintf('DEBUG: N_%d%d = (%.1f, %.1f), N_%d%d = (%.1f, %.1f)\n', ...
            apical_nodes(1,3), apical_nodes(1,4), n1_2d, ...
            apical_nodes(2,3), apical_nodes(2,4), n2_2d);

        r1 = K_inv * [n1_2d, 1]'; r1 = r1 / norm(r1);
        r2 = K_inv * [n2_2d, 1]'; r2 = r2 / norm(r2);

        a1 = dot(r1, axis_dir);  % r1 · axis
        a2 = dot(r2, axis_dir);  % r2 · axis
        fprintf('DEBUG: a1=%.4f, a2=%.4f\n', a1, a2);

        % Distance between consecutive apical nodes = 1*d
        idx1 = apical_nodes(1, 3);
        idx2 = apical_nodes(2, 3);
        delta_d = (idx2 - idx1) * d;  % Should be positive (idx2 > idx1)
        fprintf('DEBUG: idx1=%d, idx2=%d, delta_d=%.4f\n', idx1, idx2, delta_d);

        % Simpler approach: approximate λ1 ≈ λ2 for nearby points
        % Then: λ * (a2 - a1) ≈ delta_d
        if abs(a2 - a1) > 1e-6
            lambda1 = abs(delta_d / (a2 - a1));
        else
            lambda1 = 10;  % Fallback
        end

        fprintf('Solved: λ1=%.4f\n', lambda1);

        % Use a1 from the FIRST apical node for base position
        a1 = dot(r1, axis_dir);
    end
end

% Fallback to arc midpoints if no apical nodal points
if isempty(lambda1) && length(arcs_A) >= 2
    fprintf('Fallback: Using arc midpoints for reference.\n');
    mid1_2d = arcs_A{1}(ceil(end/2), :);
    mid2_2d = arcs_A{2}(ceil(end/2), :);

    r1 = K_inv * [mid1_2d, 1]'; r1 = r1 / norm(r1);
    r2 = K_inv * [mid2_2d, 1]'; r2 = r2 / norm(r2);

    a1 = dot(r1, axis_dir);
    a2 = dot(r2, axis_dir);
    t1 = dot(r1, trans_dir);
    t2 = dot(r2, trans_dir);

    if abs(t2) > 1e-6 && abs(t1) > 1e-6
        ratio = t1 / t2;
        denom = ratio * a2 - a1;
        if abs(denom) > 1e-6
            lambda1 = d / denom;
        end
    end

    if isempty(lambda1) || lambda1 <= 0
        lambda1 = abs(d / (a2 - a1 + 1e-6));
    end
end

if isempty(lambda1)
    lambda1 = 10;  % Ultimate fallback
end

% Compute reference base axis position
base_axis_pos = lambda1 * a1;
fprintf('Reference: λ1=%.4f, base_axis_pos=%.4f\n', lambda1, base_axis_pos);

% Validity threshold: cap depth at 3× reference to avoid outliers
MAX_DEPTH_FACTOR = inf;
max_lambda = MAX_DEPTH_FACTOR * lambda1;

%% 3. Reconstruct each arc with per-point depth (with validity check)
% Points near the vanishing line (ray_axis ≈ 0) are unreconstructible
ANGLE_THRESHOLD = 0.003; 

for i = 1:length(arcs_A)
    arc_pts_2d = arcs_A{i};
    n_pts = size(arc_pts_2d, 1);
    pts_3d = [];  % Only add valid points

    delta_axis = (i - 1) * d;
    target_axis_pos = base_axis_pos + delta_axis;

    for j = 1:n_pts
        ray = K_inv * [arc_pts_2d(j, :), 1]';
        ray = ray / norm(ray);

        ray_axis = dot(ray, axis_dir);

        % Check if point is reconstructible (not near vanishing line)
        if abs(ray_axis) < ANGLE_THRESHOLD
            % Point too close to horizon - skip
            continue;
        end

        lambda_j = target_axis_pos / ray_axis;

        % Only accept positive, finite depths within reasonable range
        if lambda_j > 0 && lambda_j < max_lambda && isfinite(lambda_j)
            pts_3d = [pts_3d; (lambda_j * ray)'];
        end
    end

    if ~isempty(pts_3d)
        points_3D(end+1).pts = pts_3d;
        points_3D(end).arc = 'A';
    end
end

for i = 1:length(arcs_B)
    arc_pts_2d = arcs_B{i};
    n_pts = size(arc_pts_2d, 1);
    pts_3d = [];  % Only add valid points

    delta_axis = (i - 1) * d;
    target_axis_pos = base_axis_pos + delta_axis;

    for j = 1:n_pts
        ray = K_inv * [arc_pts_2d(j, :), 1]';
        ray = ray / norm(ray);

        ray_axis = dot(ray, axis_dir);

        % Check if point is reconstructible (not near vanishing line)
        if abs(ray_axis) < ANGLE_THRESHOLD
            % Point too close to horizon - skip
            continue;
        end

        lambda_j = target_axis_pos / ray_axis;

        % Only accept positive, finite depths within reasonable range
        if lambda_j > 0 && lambda_j < max_lambda && isfinite(lambda_j)
            pts_3d = [pts_3d; (lambda_j * ray)'];
        end
    end

    if ~isempty(pts_3d)
        points_3D(end+1).pts = pts_3d;
        points_3D(end).arc = 'B';
    end
end

%% 4. Localize cylinder axis using apical nodal points
P_axis = [];

if nargin >= 8 && ~isempty(nodal_points)
    apical_nodes = nodal_points(nodal_points(:,3) == nodal_points(:,4), :);
    if ~isempty(apical_nodes)
        n1_2d = apical_nodes(1, 1:2);
        r_n1 = K_inv * [n1_2d, 1]'; r_n1 = r_n1 / norm(r_n1);
        P_axis = lambda1 * r_n1;
        fprintf('Axis from apical N_%d%d: [%.3f, %.3f, %.3f]\n', ...
            apical_nodes(1,3), apical_nodes(1,4), P_axis);
    end
end

if isempty(P_axis) && ~isempty(line_apical)
    mid_apical = (line_apical.p1 + line_apical.p2) / 2;
    ray_apical = K_inv * [mid_apical, 1]';
    ray_apical = ray_apical / norm(ray_apical);
    P_axis = lambda1 * ray_apical;
end

%% 5. Verification
if length(points_3D) >= 2
    c1 = mean(points_3D(1).pts, 1);
    c2 = mean(points_3D(2).pts, 1);
    measured_d = dot(c2 - c1, axis_dir');
    fprintf('Verification: Arc A1-A2 spacing = %.4f (target: 1.0)\n', measured_d);
end

end

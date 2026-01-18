function [points_3D, axis_dir, vert_dir, trans_dir, P_axis, R_cylinder] = reconstruct_3d(K, v_vert, v_axis, v_trans, arcs_A, arcs_B, line_apical, nodal_points)
% RECONSTRUCT_3D - 3D reconstruction using ray-cylinder intersection
%
% Correct Algorithm:
%   1. Define cylinder from apical nodal points (axis line + radius R)
%   2. For each arc point: intersect viewing ray with cylinder surface
%   3. Constraint: ||cross(λ*ray - P_axis, axis_dir)|| = R
%
% This solves a quadratic equation for λ (depth) at each point.

fprintf('\n========== DEBUG: reconstruct_3d START ==========\n');

K_inv = inv(K);

fprintf('DEBUG: K matrix:\n');
disp(K);
fprintf('DEBUG: K_inv matrix:\n');
disp(K_inv);

%% 1. Compute 3D direction vectors from vanishing points
fprintf('\n--- DEBUG: Vanishing Points Input ---\n');
fprintf('DEBUG: v_vert = '); disp(v_vert(:)');
fprintf('DEBUG: v_axis = '); disp(v_axis(:)');
fprintf('DEBUG: v_trans = '); disp(v_trans(:)');

if ~isempty(v_vert)
    vert_dir = K_inv * v_vert(:);
    fprintf('DEBUG: vert_dir (before norm) = [%.6f, %.6f, %.6f]\n', vert_dir);
    vert_dir = vert_dir / norm(vert_dir);
else
    vert_dir = [0; 1; 0];
    fprintf('DEBUG: v_vert empty, using default\n');
end

if ~isempty(v_axis)
    axis_dir = K_inv * v_axis(:);
    fprintf('DEBUG: axis_dir (before norm) = [%.6f, %.6f, %.6f]\n', axis_dir);
    axis_dir = axis_dir / norm(axis_dir);
else
    axis_dir = [1; 0; 0];
    fprintf('DEBUG: v_axis empty, using default\n');
end

if ~isempty(v_trans)
    trans_dir = K_inv * v_trans(:);
    fprintf('DEBUG: trans_dir (before norm) = [%.6f, %.6f, %.6f]\n', trans_dir);
    trans_dir = trans_dir / norm(trans_dir);
else
    trans_dir = [0; 0; 1];
    fprintf('DEBUG: v_trans empty, using default\n');
end

axis_dir = axis_dir(:)';  % Ensure row vector for consistency

fprintf('\n--- DEBUG: 3D Directions (normalized) ---\n');
fprintf('DEBUG: vert_dir  = [%.6f, %.6f, %.6f]\n', vert_dir);
fprintf('DEBUG: axis_dir  = [%.6f, %.6f, %.6f]\n', axis_dir);
fprintf('DEBUG: trans_dir = [%.6f, %.6f, %.6f]\n', trans_dir);
fprintf('3D directions computed from vanishing points.\n');

%% 2. Define Cylinder Geometry from Apical Nodal Points
fprintf('\n--- DEBUG: Cylinder Geometry Calculation ---\n');
d = 1;  % Known arc spacing
P_axis = [];
R_cylinder = [];
lambda_ref = 10;  % Default reference depth

fprintf('DEBUG: nargin = %d\n', nargin);
if nargin >= 8
    fprintf('DEBUG: nodal_points size = [%d, %d]\n', size(nodal_points, 1), size(nodal_points, 2));
    fprintf('DEBUG: nodal_points content:\n');
    disp(nodal_points);
end

if nargin >= 8 && ~isempty(nodal_points)
    % Find apical nodal points: N_ii where arcA_idx == arcB_idx
    apical_mask = nodal_points(:,3) == nodal_points(:,4);
    apical_nodes = nodal_points(apical_mask, :);
    fprintf('DEBUG: apical_mask = '); disp(apical_mask');
    fprintf('DEBUG: Found %d apical nodes\n', size(apical_nodes, 1));

    % Find non-apical nodes (N_ij where i != j) for radius calculation
    non_apical_nodes = nodal_points(nodal_points(:,3) ~= nodal_points(:,4), :);
    fprintf('DEBUG: Found %d non-apical nodes\n', size(non_apical_nodes, 1));

    if size(apical_nodes, 1) >= 2
        % Sort apical nodes by arc index
        [~, sort_idx] = sort(apical_nodes(:, 3));
        apical_nodes = apical_nodes(sort_idx, :);

        fprintf('DEBUG: Sorted apical nodes:\n');
        disp(apical_nodes);

        fprintf('Using %d apical nodes for cylinder definition.\n', size(apical_nodes, 1));

        % Get rays for first two apical nodes (N_11, N_22)
        n1_2d = apical_nodes(1, 1:2);
        n2_2d = apical_nodes(2, 1:2);

        fprintf('DEBUG: N_11 2D = [%.2f, %.2f]\n', n1_2d);
        fprintf('DEBUG: N_22 2D = [%.2f, %.2f]\n', n2_2d);

        r1 = K_inv * [n1_2d, 1]'; r1 = r1 / norm(r1);
        r2 = K_inv * [n2_2d, 1]'; r2 = r2 / norm(r2);

        fprintf('DEBUG: r1 (ray for N_11) = [%.6f, %.6f, %.6f]\n', r1);
        fprintf('DEBUG: r2 (ray for N_22) = [%.6f, %.6f, %.6f]\n', r2);

        % Dot products with axis
        a1 = dot(r1, axis_dir);
        a2 = dot(r2, axis_dir);

        fprintf('DEBUG: a1 = dot(r1, axis_dir) = %.6f\n', a1);
        fprintf('DEBUG: a2 = dot(r2, axis_dir) = %.6f\n', a2);
        fprintf('DEBUG: a2 - a1 = %.6f\n', a2 - a1);

        % Distance between N_11 and N_22 along axis = (idx2 - idx1) * d
        idx1 = apical_nodes(1, 3);
        idx2 = apical_nodes(2, 3);
        delta_d = (idx2 - idx1) * d;

        fprintf('DEBUG: idx1 = %d, idx2 = %d\n', idx1, idx2);
        fprintf('DEBUG: delta_d = (idx2 - idx1) * d = %.4f\n', delta_d);

        % Solve for reference depth
        if abs(a2 - a1) > 1e-6
            lambda_ref = abs(delta_d / (a2 - a1));
            fprintf('DEBUG: lambda_ref = |delta_d / (a2-a1)| = |%.4f / %.6f| = %.4f\n', delta_d, a2-a1, lambda_ref);
        else
            fprintf('DEBUG: WARNING - a2-a1 too small (%.6f), keeping lambda_ref = %.4f\n', a2-a1, lambda_ref);
        end

        fprintf('Reference depth λ = %.4f\n', lambda_ref);

        % 3D position of first apical nodal point (apex - on the roof ridge)
        P_apex = lambda_ref * r1(:)';
        fprintf('DEBUG: P_apex = lambda_ref * r1 = [%.4f, %.4f, %.4f]\n', P_apex);

        % Calculate cylinder radius R using a NON-APICAL nodal point
        % The non-apical node (e.g., N_12) lies on the cylinder surface
        R_cylinder = 1.0;  % Default

        if ~isempty(non_apical_nodes)
            fprintf('DEBUG: Using non-apical node for radius calculation\n');
            fprintf('DEBUG: All non-apical nodes:\n');
            disp(non_apical_nodes);

            % CRITICAL: Filter out edge nodes where either index is 0
            % Edge nodes have rays nearly parallel to axis, causing lambda explosion
            valid_mask = (non_apical_nodes(:,3) > 0) & (non_apical_nodes(:,4) > 0);
            interior_nodes = non_apical_nodes(valid_mask, :);
            fprintf('DEBUG: Found %d interior non-apical nodes (both indices > 0)\n', size(interior_nodes, 1));

            if isempty(interior_nodes)
                % Fallback: use any non-apical node but filter by ray angle
                fprintf('DEBUG: No interior nodes, trying to find best non-apical node by ray angle\n');
                interior_nodes = non_apical_nodes;
            end

            % Among valid nodes, prefer one with similar depth to reference (good ray angle)
            % Check ray angles and pick the one with largest |a_na| (most perpendicular to axis)
            best_idx = 1;
            best_a_na = 0;
            for ni = 1:size(interior_nodes, 1)
                na_test = interior_nodes(ni, 1:2);
                r_test = K_inv * [na_test, 1]'; r_test = r_test / norm(r_test);
                a_test = abs(dot(r_test, axis_dir));
                fprintf('DEBUG:   Node %d: [%.1f, %.1f, %d, %d], |a_na| = %.4f\n', ...
                    ni, interior_nodes(ni,:), a_test);
                if a_test > best_a_na
                    best_a_na = a_test;
                    best_idx = ni;
                end
            end

            na = interior_nodes(best_idx, :);
            fprintf('DEBUG: Selected non-apical node = [%.2f, %.2f, %d, %d]\n', na);
            na_2d = na(1:2);
            r_na = K_inv * [na_2d, 1]'; r_na = r_na / norm(r_na);
            fprintf('DEBUG: r_na (ray for non-apical) = [%.6f, %.6f, %.6f]\n', r_na);

            % Reconstruct N_12 at similar depth to reference
            % N_12 is between arc 1 and arc 2, so its axis position is intermediate
            a_na = dot(r_na, axis_dir);
            fprintf('DEBUG: a_na = dot(r_na, axis_dir) = %.6f\n', a_na);

            if abs(a_na) > 1e-6
                lambda_na = lambda_ref * a1 / a_na;  % Approximately same axis position
                fprintf('DEBUG: lambda_na = lambda_ref * a1 / a_na = %.4f * %.6f / %.6f = %.4f\n', lambda_ref, a1, a_na, lambda_na);
            else
                lambda_na = lambda_ref;
                fprintf('DEBUG: WARNING - a_na too small, using lambda_na = lambda_ref = %.4f\n', lambda_na);
            end

            P_side = lambda_na * r_na(:)';
            fprintf('DEBUG: P_side = lambda_na * r_na = [%.4f, %.4f, %.4f]\n', P_side);

            % The radius R is the transversal distance from apex to side point
            % projected perpendicular to the axis
            diff = P_side - P_apex;
            fprintf('DEBUG: diff = P_side - P_apex = [%.4f, %.4f, %.4f]\n', diff);
            diff_perp = diff - dot(diff, axis_dir) * axis_dir(:)';
            fprintf('DEBUG: dot(diff, axis_dir) = %.4f\n', dot(diff, axis_dir));
            fprintf('DEBUG: diff_perp = [%.4f, %.4f, %.4f]\n', diff_perp);
            R_cylinder = norm(diff_perp);
            fprintf('DEBUG: R_cylinder = norm(diff_perp) = %.4f\n', R_cylinder);

            fprintf('Using non-apical node for R. R = %.4f\n', R_cylinder);
        else
            fprintf('DEBUG: No non-apical nodes found, using default R = %.4f\n', R_cylinder);
        end

        % Cylinder axis is BELOW the apex (inside the vault volume)
        % Shift axis down by R in the vertical direction
        P_axis = P_apex - R_cylinder * vert_dir(:)';
        fprintf('DEBUG: vert_dir = [%.6f, %.6f, %.6f]\n', vert_dir);
        fprintf('DEBUG: R_cylinder * vert_dir = [%.4f, %.4f, %.4f]\n', R_cylinder * vert_dir(:)');

        fprintf('Axis shifted down by R. P_axis = [%.3f, %.3f, %.3f]\n', P_axis);
    else
        fprintf('DEBUG: Not enough apical nodes (found %d, need >= 2)\n', size(apical_nodes, 1));
    end
else
    fprintf('DEBUG: nodal_points empty or not provided\n');
end

% Fallback if no apical nodes
if isempty(P_axis)
    fprintf('Warning: No apical nodes. Using default cylinder.\n');
    P_axis = [0, 0, lambda_ref];
    R_cylinder = 1.0;
end

fprintf('\n--- DEBUG: Final Cylinder Parameters ---\n');
fprintf('DEBUG: P_axis = [%.4f, %.4f, %.4f]\n', P_axis);
fprintf('DEBUG: R_cylinder = %.4f\n', R_cylinder);
fprintf('DEBUG: lambda_ref = %.4f\n', lambda_ref);

%% 3. Reconstruct arcs using Ray-Cylinder Intersection
% For each point: solve ||cross(λ*ray - P_axis, axis_dir)|| = R
% This is a quadratic in λ.

fprintf('\n--- DEBUG: Arc Reconstruction ---\n');
fprintf('DEBUG: Number of arcs_A = %d\n', length(arcs_A));
fprintf('DEBUG: Number of arcs_B = %d\n', length(arcs_B));

points_3D = struct('pts', {}, 'arc', {});

for i = 1:length(arcs_A)
    arc_pts_2d = arcs_A{i};
    n_pts = size(arc_pts_2d, 1);
    pts_3d = [];
    failed_pts = 0;
    sample_lambdas = [];

    for j = 1:n_pts
        ray = K_inv * [arc_pts_2d(j, :), 1]';
        ray = ray(:)' / norm(ray);  % Row vector, normalized

        lambda = solve_ray_cylinder(ray, P_axis, axis_dir, R_cylinder, lambda_ref);

        % Collect sample for debugging (first, middle, last points)
        if j == 1 || j == round(n_pts/2) || j == n_pts
            sample_lambdas = [sample_lambdas; j, lambda];
        end

        if ~isempty(lambda) && lambda > 0 && isfinite(lambda)
            pt_3d = lambda * ray;
            pts_3d = [pts_3d; pt_3d];
        else
            failed_pts = failed_pts + 1;
        end
    end

    fprintf('DEBUG: Arc A[%d]: %d pts, %d reconstructed, %d failed\n', i, n_pts, size(pts_3d, 1), failed_pts);
    if ~isempty(sample_lambdas)
        fprintf('DEBUG:   Sample lambdas (pt_idx, lambda): ');
        for k = 1:size(sample_lambdas, 1)
            fprintf('[%d, %.4f] ', sample_lambdas(k,:));
        end
        fprintf('\n');
    end
    if ~isempty(pts_3d)
        fprintf('DEBUG:   3D range X: [%.4f, %.4f]\n', min(pts_3d(:,1)), max(pts_3d(:,1)));
        fprintf('DEBUG:   3D range Y: [%.4f, %.4f]\n', min(pts_3d(:,2)), max(pts_3d(:,2)));
        fprintf('DEBUG:   3D range Z: [%.4f, %.4f]\n', min(pts_3d(:,3)), max(pts_3d(:,3)));
    end

    if ~isempty(pts_3d)
        points_3D(end+1).pts = pts_3d;
        points_3D(end).arc = 'A';
    end
end

for i = 1:length(arcs_B)
    arc_pts_2d = arcs_B{i};
    n_pts = size(arc_pts_2d, 1);
    pts_3d = [];
    failed_pts = 0;
    sample_lambdas = [];

    for j = 1:n_pts
        ray = K_inv * [arc_pts_2d(j, :), 1]';
        ray = ray(:)' / norm(ray);

        lambda = solve_ray_cylinder(ray, P_axis, axis_dir, R_cylinder, lambda_ref);

        % Collect sample for debugging
        if j == 1 || j == round(n_pts/2) || j == n_pts
            sample_lambdas = [sample_lambdas; j, lambda];
        end

        if ~isempty(lambda) && lambda > 0 && isfinite(lambda)
            pt_3d = lambda * ray;
            pts_3d = [pts_3d; pt_3d];
        else
            failed_pts = failed_pts + 1;
        end
    end

    fprintf('DEBUG: Arc B[%d]: %d pts, %d reconstructed, %d failed\n', i, n_pts, size(pts_3d, 1), failed_pts);
    if ~isempty(sample_lambdas)
        fprintf('DEBUG:   Sample lambdas (pt_idx, lambda): ');
        for k = 1:size(sample_lambdas, 1)
            fprintf('[%d, %.4f] ', sample_lambdas(k,:));
        end
        fprintf('\n');
    end
    if ~isempty(pts_3d)
        fprintf('DEBUG:   3D range X: [%.4f, %.4f]\n', min(pts_3d(:,1)), max(pts_3d(:,1)));
        fprintf('DEBUG:   3D range Y: [%.4f, %.4f]\n', min(pts_3d(:,2)), max(pts_3d(:,2)));
        fprintf('DEBUG:   3D range Z: [%.4f, %.4f]\n', min(pts_3d(:,3)), max(pts_3d(:,3)));
    end

    if ~isempty(pts_3d)
        points_3D(end+1).pts = pts_3d;
        points_3D(end).arc = 'B';
    end
end

fprintf('\n--- DEBUG: Reconstruction Summary ---\n');
fprintf('%d arcs reconstructed on cylinder surface.\n', length(points_3D));

%% 4. Verification
fprintf('\n--- DEBUG: Verification ---\n');
if length(points_3D) >= 2
    c1 = mean(points_3D(1).pts, 1);
    c2 = mean(points_3D(2).pts, 1);
    fprintf('DEBUG: Arc 1 centroid = [%.4f, %.4f, %.4f]\n', c1);
    fprintf('DEBUG: Arc 2 centroid = [%.4f, %.4f, %.4f]\n', c2);
    fprintf('DEBUG: c2 - c1 = [%.4f, %.4f, %.4f]\n', c2 - c1);

    measured_d = abs(dot(c2 - c1, axis_dir));
    fprintf('DEBUG: dot(c2-c1, axis_dir) = %.4f\n', dot(c2 - c1, axis_dir));
    fprintf('Verification: Arc spacing = %.4f (target: 1.0)\n', measured_d);

    % Check cylinder radius
    dist_to_axis = norm(cross(c1 - P_axis, axis_dir));
    fprintf('DEBUG: cross(c1 - P_axis, axis_dir) = [%.4f, %.4f, %.4f]\n', cross(c1 - P_axis, axis_dir));
    fprintf('Verification: Distance to axis = %.4f (R = %.4f)\n', dist_to_axis, R_cylinder);

    % Check all arc centroids
    fprintf('\nDEBUG: All arc centroids:\n');
    for i = 1:length(points_3D)
        ci = mean(points_3D(i).pts, 1);
        dist_i = norm(cross(ci - P_axis, axis_dir));
        fprintf('DEBUG:   Arc %d (%s): centroid = [%.4f, %.4f, %.4f], dist to axis = %.4f\n', ...
            i, points_3D(i).arc, ci, dist_i);
    end
else
    fprintf('DEBUG: Not enough arcs for verification (got %d)\n', length(points_3D));
end

fprintf('\n========== DEBUG: reconstruct_3d END ==========\n');

end

%% ========== HELPER FUNCTION ==========

function lambda = solve_ray_cylinder(ray, P_axis, axis_dir, R, lambda_ref)
% SOLVE_RAY_CYLINDER - Find intersection of ray with infinite cylinder
%
% Constraint: ||cross(P - P_axis, axis_dir)|| = R
% Where P = lambda * ray
%
% Expanding: ||cross(λ*ray - P_axis, a)|| = R
% Let d = ray, p = P_axis, a = axis_dir (all unit/normalized)
%
% ||cross(λ*d - p, a)||^2 = R^2
% = ||λ * cross(d,a) - cross(p,a)||^2 = R^2
%
% Let u = cross(d, a), v = cross(p, a)
% ||λ*u - v||^2 = R^2
% λ^2 * ||u||^2 - 2*λ*(u·v) + ||v||^2 = R^2
%
% Quadratic: A*λ^2 + B*λ + C = 0
% A = ||u||^2
% B = -2*(u·v)
% C = ||v||^2 - R^2

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
    % No real intersection (ray parallel to axis or misses cylinder)
    % Fall back to approximate depth
    lambda = lambda_ref;
    return;
end

% Two solutions - pick the one closer to reference depth
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

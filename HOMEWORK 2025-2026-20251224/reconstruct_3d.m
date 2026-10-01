function [points_3D, axis_dir, vert_dir, trans_dir, P_axis, R_cylinder] = reconstruct_3d(K, v_vert, v_axis, v_trans, arcs_A, arcs_B, line_apical, nodal_points)
% RECONSTRUCT_3D - 3D reconstruction using ray-cylinder intersection
%
% Algorithm:
%   1. Define cylinder from apical nodal points (axis line + radius R)
%   2. For each arc point: intersect viewing ray with cylinder surface
%   3. Constraint: ||cross(λ*ray - P_axis, axis_dir)|| = R
%
% This solves a quadratic equation for λ (depth) at each point.

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

axis_dir = axis_dir(:)';  % Ensure row vector for consistency

%% 2. Define Cylinder Geometry from Apical Nodal Points
d = 1;  % Known arc spacing (metric constraint)
P_axis = [];
R_cylinder = [];
lambda_ref = [];

if nargin >= 8 && ~isempty(nodal_points)
    % Find apical nodal points: N_ii where arcA_idx == arcB_idx
    apical_mask = nodal_points(:,3) == nodal_points(:,4);
    apical_nodes = nodal_points(apical_mask, :);

    % Find non-apical nodes (N_ij where i != j) for radius calculation
    non_apical_nodes = nodal_points(~apical_mask, :);

    if size(apical_nodes, 1) >= 2
        % Sort apical nodes by arc index
        [~, sort_idx] = sort(apical_nodes(:, 3));
        apical_nodes = apical_nodes(sort_idx, :);

        % Get rays for first two apical nodes (N_11, N_22)
        n1_2d = apical_nodes(1, 1:2);
        n2_2d = apical_nodes(2, 1:2);

        r1 = K_inv * [n1_2d, 1]'; r1 = r1 / norm(r1);
        r2 = K_inv * [n2_2d, 1]'; r2 = r2 / norm(r2);

        % Dot products with axis
        a1 = dot(r1, axis_dir);
        a2 = dot(r2, axis_dir);

        % Distance between N_11 and N_22 along axis = (idx2 - idx1) * d
        idx1 = apical_nodes(1, 3);
        idx2 = apical_nodes(2, 3);
        delta_d = (idx2 - idx1) * d;

        % Compute reference depth from the d=1 constraint
        if abs(a2 - a1) > 1e-6
            lambda_ref = abs(delta_d / (a2 - a1));
        else
            lambda_ref = 5.0;  % Fallback
        end

        % 3D position of first apical nodal point (apex - on the roof ridge)
        P_apex = lambda_ref * r1(:)';

        % Calculate cylinder radius R using ALL non-apical nodes (robust median)
        R_cylinder = d;  % Default

        if ~isempty(non_apical_nodes)
            n_nodes = size(non_apical_nodes, 1);
            lambdas = zeros(n_nodes, 1);

            % Compute lambda for all non-apical nodes
            for ni = 1:n_nodes
                na_test = non_apical_nodes(ni, 1:2);
                r_test = K_inv * [na_test, 1]'; r_test = r_test / norm(r_test);
                a_test = dot(r_test, axis_dir);
                if abs(a_test) > 1e-6
                    lambdas(ni) = lambda_ref * a1 / a_test;
                else
                    lambdas(ni) = Inf;
                end
            end

            % Outlier detection using median absolute deviation (MAD)
            finite_lambdas = lambdas(isfinite(lambdas));
            if ~isempty(finite_lambdas)
                med_lambda = median(finite_lambdas);
                mad_lambda = median(abs(finite_lambdas - med_lambda));
                threshold = max(3 * mad_lambda, 0.5 * lambda_ref);
                inlier_mask = abs(lambdas - med_lambda) < threshold & isfinite(lambdas);
            else
                inlier_mask = true(n_nodes, 1);
            end

            % Compute R from all inlier nodes, take median
            R_estimates = [];
            for ni = 1:n_nodes
                if inlier_mask(ni)
                    na_i = non_apical_nodes(ni, 1:2);
                    r_i = K_inv * [na_i, 1]'; r_i = r_i / norm(r_i);
                    P_side_i = lambdas(ni) * r_i(:)';
                    diff_i = P_side_i - P_apex;
                    diff_perp_i = diff_i - dot(diff_i, axis_dir) * axis_dir(:)';
                    R_estimates = [R_estimates; norm(diff_perp_i)];
                end
            end

            if ~isempty(R_estimates)
                R_cylinder = median(R_estimates);
            end
        end

        % Cylinder axis is BELOW the apex (inside the vault volume)
        P_axis = P_apex - R_cylinder * vert_dir(:)';
    end
end

% Fallback if no apical nodes
if isempty(P_axis)
    warning('No apical nodes found. Using default cylinder.');
    if isempty(lambda_ref), lambda_ref = 5.0; end
    P_axis = [0, 0, lambda_ref];
    R_cylinder = 1.0;
end

%% 3. Reconstruct arcs using Ray-Cylinder Intersection
points_3D = struct('pts', {}, 'arc', {});

for i = 1:length(arcs_A)
    arc_pts_2d = arcs_A{i};
    n_pts = size(arc_pts_2d, 1);
    pts_3d = [];

    for j = 1:n_pts
        ray = K_inv * [arc_pts_2d(j, :), 1]';
        ray = ray(:)' / norm(ray);
        lambda = solve_ray_cylinder(ray, P_axis, axis_dir, R_cylinder, lambda_ref);
        if ~isempty(lambda) && lambda > 0 && isfinite(lambda)
            pts_3d = [pts_3d; lambda * ray];
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
    pts_3d = [];

    for j = 1:n_pts
        ray = K_inv * [arc_pts_2d(j, :), 1]';
        ray = ray(:)' / norm(ray);
        lambda = solve_ray_cylinder(ray, P_axis, axis_dir, R_cylinder, lambda_ref);
        if ~isempty(lambda) && lambda > 0 && isfinite(lambda)
            pts_3d = [pts_3d; lambda * ray];
        end
    end

    if ~isempty(pts_3d)
        points_3D(end+1).pts = pts_3d;
        points_3D(end).arc = 'B';
    end
end

%% 4. Verification (brief output)
if length(points_3D) >= 2
    c1 = mean(points_3D(1).pts, 1);
    c2 = mean(points_3D(2).pts, 1);
    measured_d = abs(dot(c2 - c1, axis_dir));
    fprintf('Arc spacing: %.3f (target: 1.0), R: %.3f\n', measured_d, R_cylinder);
end

end

%% ========== HELPER FUNCTION ==========

function lambda = solve_ray_cylinder(ray, P_axis, axis_dir, R, lambda_ref)
% SOLVE_RAY_CYLINDER - Find intersection of ray with infinite cylinder
%
% Solves: ||cross(λ*ray - P_axis, axis_dir)|| = R
% This is a quadratic equation in λ.

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

candidates = [lambda1, lambda2];
candidates = candidates(candidates > 0);

if isempty(candidates)
    lambda = lambda_ref;
else
    [~, idx] = min(abs(candidates - lambda_ref));
    lambda = candidates(idx);
end

end

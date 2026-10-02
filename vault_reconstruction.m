% vault_reconstruction.m - Image Analysis and Computer Vision Homework 2025-2026
% Student: Karim Negm
% Project: 3D Geometric Reconstruction of a Cylindrical Vault
% Description: This script performs a complete pipeline for metric 3D reconstruction
%              of the San Maurizio church vault from a single uncalibrated image.
%              Main stages: vanishing point estimation, stratified rectification,
%              camera calibration, and 3D reconstruction via cylindrical constraints.
%
% Revised after submission (October 2026): the aspect ratio of the metric rectification
% now comes from cross-section circles built from mirrored arc pairs, the nodal points
% are intersected in the image, and the cylinder is fitted to those cross-sections.

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

%% Part 3b: Stratified Metric Rectification (Circular Cross-Section)
% Two constraints: (1) Orthogonality fixes skew, (2) Circular profile fixes aspect ratio.
fprintf('\n--- Part 3b: Metric Rectification (Circular Cross-Section) ---\n');

% CONSTRAINT 1: Orthogonality (fixes skew)
H_ortho = [S_metric_init, [0;0]; 0 0 1];
H_affine_ortho = H_ortho * H_affine;

% CONSTRAINT 2: Circular cross-section (fixes aspect ratio λ)
% Arcs a_i and b_j are mirror images in the plane pi_ij perpendicular to the axis.
% A line through the axis vanishing point is the image of a line parallel to the axis.
% Where it meets a_i (at p) and b_j (at q), the two 3D points are mirror images, so
% their 3D midpoint lies in pi_ij and on the cylinder: it is a point of the circular
% cross-section. Its image is the harmonic conjugate of the axis vanishing point with
% respect to p and q. Each pair of arcs therefore gives the image of one cross-section
% circle. The rectifying homography is the same for all planes perpendicular to the
% axis, so all of these circles must become circles under the same λ.
% Arcs that cross the vanishing line are left out: their planes pass close to the camera.
usable_A = find(cellfun(@(a) ~crosses_line(a, vanishing_line), arcs_A));
usable_B = find(cellfun(@(b) ~crosses_line(b, vanishing_line), arcs_B));

sections = struct('i', {}, 'j', {}, 'pts', {});
for i = usable_A(:)'
    for j = usable_B(:)'
        mid = mirror_midpoints(smooth_arc(arcs_A{i}), smooth_arc(arcs_B{j}), v_axis_pixel(1:2));
        if size(mid, 1) >= 8
            sections(end+1) = struct('i', i, 'j', j, 'pts', mid);
        end
    end
end
n_sec = numel(sections);
fprintf('Cross-section circles recovered from %d mirrored arc pairs.\n', n_sec);

% In the affine-orthogonal frame (x vertical, y transversal) the circle of cross-section k is
%   x^2 + λ^2 y^2 + a_k x + b_k y + c_k = 0,
% with the same λ for every k. This is linear in the unknowns and solved for all k at once.
best_lambda = fit_aspect_ratio(sections, H_affine_ortho);
lambda_estimates = arrayfun(@(s) fit_aspect_ratio(s, H_affine_ortho), sections);
if isnan(best_lambda)
    error('The cross-section points do not define an ellipse.');
end
fprintf('λ per arc pair: range=[%.3f, %.3f]\n', min(lambda_estimates), max(lambda_estimates));
fprintf('Constraint 2 (Circular): λ = %.4f (joint fit over %d cross-sections)\n', best_lambda, n_sec);

% Construct final H_metric
H_scale = [1, 0, 0; 0, best_lambda, 0; 0, 0, 1];
H_metric = H_scale * H_ortho;
fprintf('Metric rectification complete: AR=%.4f\n', best_lambda);

% Final combined homography
H_rectify = H_metric * H_affine;

% --- Auto-Rotation Correction ---
% The metric rectification may produce an arbitrary orientation.
% We rotate so that the vertical vanishing point direction becomes truly vertical.
v_vert_rect = H_rectify * v_vert_pixel(:);
v_vert_rect = v_vert_rect / v_vert_rect(3);  % Normalize

% Compute the direction in rectified space (as a direction at infinity, use first 2 components)
v_vert_affine = H_rectify * v_vert_pixel(:);
vert_dir = v_vert_affine(1:2);  % Direction vector
vert_dir = vert_dir / norm(vert_dir);

% We want vertical to point "up" (i.e., along -Y in image coordinates, or [0, -1])
target_dir = [0; -1];

% Compute rotation angle to align vert_dir with target_dir
angle = atan2(vert_dir(1)*target_dir(2) - vert_dir(2)*target_dir(1), ...
    vert_dir(1)*target_dir(1) + vert_dir(2)*target_dir(2));

% Build rotation matrix
R_align = [cos(angle), -sin(angle), 0; ...
    sin(angle),  cos(angle), 0; ...
    0,           0,          1];

% Apply rotation to rectification homography
H_rectify = R_align * H_rectify;

fprintf('Applied rotation correction of %.1f degrees to align vertical direction.\n', rad2deg(angle));

% --- Vertical Flip Correction ---
% Check if the image is upside-down by testing where a top point goes
% Transform a point from the top of the original image and bottom
top_pt = H_rectify * [cols/2; 1; 1];  % Top center of original
bot_pt = H_rectify * [cols/2; rows; 1];  % Bottom center of original
top_pt = top_pt(1:2) / top_pt(3);
bot_pt = bot_pt(1:2) / bot_pt(3);

% If top point has larger Y than bottom point, the image is flipped
if top_pt(2) > bot_pt(2)
    % Apply vertical flip (negate Y)
    V_flip = [1, 0, 0; 0, -1, 0; 0, 0, 1];
    H_rectify = V_flip * H_rectify;
    fprintf('Applied vertical flip correction.\n');
end

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

%% Part 4: Camera Intrinsic Calibration (IAC Method)
% Estimate K using the Image of the Absolute Conic (omega).
fprintf('\n--- Computing Camera Calibration matrix K ---\n');

% The equations are written in the normalized coordinates of Part 2, so that all
% entries of the linear system have comparable size.
T_pixel_to_norm = inv(T_norm_to_pixel);

% Build constraints from vanishing point orthogonality
vps = {T_pixel_to_norm * v_vert_pixel(:), T_pixel_to_norm * v_axis_pixel(:), T_pixel_to_norm * v_trans_pixel(:)};
pairs = [1 2; 1 3; 2 3];
A_iac = [];
for k = 1:size(pairs, 1)
    A_iac = [A_iac; iac_row(vps{pairs(k,1)}, vps{pairs(k,2)})];
end

% Add constraints from rectification homography
H_world_to_img = T_pixel_to_norm * inv(H_rectify);
h1 = H_world_to_img(:,1); h2 = H_world_to_img(:,2);
A_iac = [A_iac; iac_row(h1, h2); iac_row(h1, h1) - iac_row(h2, h2)];
A_iac = A_iac ./ vecnorm(A_iac, 2, 2);

% Solve for omega and extract K
[~, ~, V_calib] = svd(A_iac);
x = V_calib(:,end);
omega = [x(1), 0, x(3); 0, x(2), x(4); x(3), x(4), x(5)];
if omega(1,1) < 0, omega = -omega; end

C = chol(omega, 'upper');
K_norm = inv(C);
K_norm = K_norm / K_norm(3,3);
K = T_norm_to_pixel * K_norm;

fprintf('K calibrated: fx=%.1f, fy=%.1f, pp=(%.1f,%.1f)\n', K(1,1), K(2,2), K(1,3), K(2,3));

% Cross-check that does not use the rectification: with square pixels, the three
% orthogonal vanishing points alone fix the principal point (the orthocentre of
% their triangle) and the focal length.
pp_sq = ([v_axis_pixel(1:2) - v_trans_pixel(1:2); v_vert_pixel(1:2) - v_trans_pixel(1:2)] \ ...
    [(v_axis_pixel(1:2) - v_trans_pixel(1:2)) * v_vert_pixel(1:2)'; ...
    (v_vert_pixel(1:2) - v_trans_pixel(1:2)) * v_axis_pixel(1:2)'])';
f_sq = sqrt(-(v_vert_pixel(1:2) - pp_sq) * (v_axis_pixel(1:2) - pp_sq)');
fprintf('Square-pixel cross-check from the vanishing points alone: f=%.1f, pp=(%.1f,%.1f)\n', ...
    f_sq, pp_sq(1), pp_sq(2));

%% Part 5: Determination of Nodal Points
% Nodal points are the critical intersection points of the vault's arcs.
% Finding these points allows us to define the skeleton of the reconstruction.
fprintf('\n--- Searching for Nodal Intersection Points ---\n');

% Find intersections between all pairs of arcs, in the image. (In the rectified plane
% the arcs that cross the vanishing line are torn apart, which creates false intersections.)
nodal_points = get_nodal_points(arcs_A, arcs_B);

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

% Get apical nodes (where A_index == B_index), sorted by index
apical_mask = nodal_points(:, 3) == nodal_points(:, 4);
apical_nodes = sortrows(nodal_points(apical_mask, :), 3);
n_ap = size(apical_nodes, 1);

% Compute the depths of the apical nodes via linear triangulation.
% The apical nodes lie on one line parallel to the axis, d apart, so for node k
%   lambda_k * ray_k - lambda_1 * ray_1 = (idx_k - idx_1) * d * axis_direction.
% All depths are solved at once in the least-squares sense.
rays_ap = zeros(3, n_ap);
for k = 1:n_ap
    ray_k = K_inv * [apical_nodes(k, 1:2), 1]';
    rays_ap(:, k) = ray_k / norm(ray_k);
end

A_tri = zeros(3 * (n_ap - 1), n_ap);
b_tri = zeros(3 * (n_ap - 1), 1);
for k = 2:n_ap
    eq = 3 * (k - 2) + (1:3);
    A_tri(eq, 1) = -rays_ap(:, 1);
    A_tri(eq, k) = rays_ap(:, k);
    b_tri(eq) = (apical_nodes(k, 3) - apical_nodes(1, 3)) * d * axis_direction(:);
end
lambdas_ap = A_tri \ b_tri;

fprintf('Triangulated depths of the apical nodes: %s\n', mat2str(lambdas_ap', 4));

% 3D position of the first apical node (apex reference).
P_apex = lambdas_ap(1) * rays_ap(:, 1)';
idx_apex = apical_nodes(1, 3);

% Cylinder radius and axis from the cross-section circles of Part 3b.
% The plane of the pair (a_i, b_j) is perpendicular to the axis at the axial coordinate
% (i + j) / 2, in units of d and with the index shift applied to j. Intersecting the
% camera rays of its midpoints with that plane gives 3D points of one cross-section.
% A circle fitted to them gives the radius and a point of the axis.
z_apex = dot(P_apex, axis_direction);
radius_estimates = zeros(n_sec, 1);
centers = zeros(n_sec, 2);
for k = 1:n_sec
    z_k = z_apex + ((sections(k).i + sections(k).j - shift) / 2 - idx_apex) * d;
    rays = K_inv * [sections(k).pts, ones(size(sections(k).pts, 1), 1)]';
    P_sec = (rays .* (z_k ./ (axis_direction * rays)))';

    % Coordinates in the cross-section plane, relative to the apex: (up, across)
    [centers(k, :), radius_estimates(k)] = fit_circle((P_sec - P_apex) * vert_direction, ...
        (P_sec - P_apex) * trans_direction);
end

R_cylinder = median(radius_estimates);
center = mean(centers, 1);
fprintf('Radius per cross-section: range=[%.3f, %.3f]\n', min(radius_estimates), max(radius_estimates));
fprintf('Estimated Cylinder Radius R: %.3f\n', R_cylinder);
fprintf('Circle center relative to the apex: %.3f up, %.3f across\n', center(1), center(2));

% Cylinder axis localization: the axis passes through the centers of the cross-sections.
P_axis = P_apex + center(1) * vert_direction(:)' + center(2) * trans_direction(:)';

% Documentation of the cylinder axis localization relative to the camera.
fprintf('\n--- Cylinder Axis Spatial Localization ---\n');
fprintf('The Axis passes through: [%.4f, %.4f, %.4f]\n', P_axis);
fprintf('Axis Direction Vector:    [%.4f, %.4f, %.4f]\n', axis_direction);
fprintf('Distance from the camera to the axis: %.3f\n', norm(P_axis - dot(P_axis, axis_direction) * axis_direction));

% Reconstruct all arc points
points_3D = {};

for i = 1:length(arcs_A)
    arc_2d = arcs_A{i};
    pts_3d = [];

    for j = 1:size(arc_2d, 1)
        ray = K_inv * [arc_2d(j, :), 1]';
        ray = ray(:)' / norm(ray);

        lambda = intersect_ray_cylinder(ray, P_axis, axis_direction, R_cylinder);

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

        lambda = intersect_ray_cylinder(ray, P_axis, axis_direction, R_cylinder);

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
fprintf('\n%d arcs reconstructed\n', length(points_3D));

% (1) Spacing of neighbouring arcs of one family, measured along the axis at equal
%     angle around the cylinder. Target: d.
spacings = [];
for family = 'AB'
    members = find(cellfun(@(s) s.arc == family, points_3D));
    for k = 1:numel(members) - 1
        [z1, th1] = cylinder_coords(points_3D{members(k)}.pts, P_axis, axis_direction, vert_direction, trans_direction);
        [z2, th2] = cylinder_coords(points_3D{members(k+1)}.pts, P_axis, axis_direction, vert_direction, trans_direction);
        th_lo = max(min(th1), min(th2));
        th_hi = min(max(th1), max(th2));
        if th_hi - th_lo < 0.1, continue; end
        th = linspace(th_lo, th_hi, 20);
        [th1, o1] = unique(th1);
        [th2, o2] = unique(th2);
        spacings(end+1) = mean(interp1(th2, z2(o2), th) - interp1(th1, z1(o1), th));
    end
end
measured_spacing = mean(spacings);
fprintf('Arc spacing verification: mean %.3f, range [%.3f, %.3f] (target: %.1f)\n', ...
    measured_spacing, min(spacings), max(spacings), d);

% (2) Reprojection of the apex line: stepping from the first apical node along the axis
%     by the node spacing must land on the other apical nodes in the image.
reproj_err = zeros(n_ap, 1);
for k = 1:n_ap
    X = P_apex + (apical_nodes(k, 3) - idx_apex) * d * axis_direction;
    x = K * X';
    reproj_err(k) = norm(x(1:2)' / x(3) - apical_nodes(k, 1:2));
end
fprintf('Apical node reprojection error: %s px\n', mat2str(reproj_err', 3));

% (3) The apex must lie on the fitted circle: its distance from the axis against R.
fprintf('Apex distance from the axis: %.3f (R = %.3f)\n', norm(center), R_cylinder);


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

% Draw vertical lines (black) with labels
for i = 1:length(lines_v)
    L = lines_v(i);
    plot([L.p1(1), L.p2(1)], [L.p1(2), L.p2(2)], 'k-', 'LineWidth', 2);
    text(mean([L.p1(1), L.p2(1)]), mean([L.p1(2), L.p2(2)]), sprintf('v%d', i), ...
        'Color', 'k', 'FontSize', 9, 'FontWeight', 'bold', 'BackgroundColor', [0.9 0.9 0.9]);
end

% Draw axis lines (green)
for i = 1:length(lines_axis)
    L = lines_axis(i);
    plot([L.p1(1), L.p2(1)], [L.p1(2), L.p2(2)], 'g-', 'LineWidth', 2);
end

% Draw transversal lines (white) with labels
for i = 1:length(lines_trans)
    L = lines_trans(i);
    plot([L.p1(1), L.p2(1)], [L.p1(2), L.p2(2)], 'w-', 'LineWidth', 2);
    text(mean([L.p1(1), L.p2(1)]), mean([L.p1(2), L.p2(2)]), sprintf('t%d', i), ...
        'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold', 'BackgroundColor', [0.2 0.2 0.2]);
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

% Plot nodal points (arc intersections)
h_nodal = plot(nodal_points(:,1), nodal_points(:,2), 'yo', 'MarkerSize', 12, ...
    'LineWidth', 2, 'MarkerFaceColor', 'y');

% Create proper legend with handles
h_vert = plot(NaN, NaN, 'k-', 'LineWidth', 2);
h_axis = plot(NaN, NaN, 'g-', 'LineWidth', 2);
h_trans = plot(NaN, NaN, 'w-', 'LineWidth', 2);
h_apic = plot(NaN, NaN, 'y-', 'LineWidth', 3);
lgd = legend([h_vert, h_axis, h_trans, h_apic, h_arcA, h_arcB, h_nodal], ...
    {'Vertical', 'Axis', 'Transversal', 'Apical', 'Arcs A', 'Arcs B', 'Nodal Points'}, ...
    'Location', 'best');
set(lgd, 'Color', [0.2 0.2 0.2], 'TextColor', 'w');

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

% Figure 7: the fitted cylinder projected back onto the photo
figure(7);
imshow(img);
title('Figure 7: Fitted Cylinder Projected onto the Image');
hold on;

z_axis_apex = dot(P_apex - P_axis, axis_direction);
theta = linspace(-115, 115, 300)' * pi / 180;
for step = -1:0.5:4
    % Cross-section of the cylinder at "step" node spacings from the first apical node
    section_center = P_axis + (z_axis_apex + step * d) * axis_direction;
    X = section_center + R_cylinder * (cos(theta) * vert_direction(:)' + sin(theta) * trans_direction(:)');
    x = project_points(K, X, cols, rows);
    if mod(step, 1) == 0
        plot(x(:,1), x(:,2), 'w-', 'LineWidth', 1.5);
    else
        plot(x(:,1), x(:,2), 'w-', 'LineWidth', 0.5);
    end
end

% Top line of the cylinder
t = linspace(-1.5, 5, 200)';
X = P_axis + (z_axis_apex + t * d) * axis_direction + R_cylinder * vert_direction(:)';
x = project_points(K, X, cols, rows);
h_top = plot(x(:,1), x(:,2), 'y-', 'LineWidth', 1.5);

h_nodes = plot(nodal_points(:,1), nodal_points(:,2), 'yo', 'MarkerSize', 10, 'LineWidth', 1.5);
h_sec = plot(NaN, NaN, 'w-', 'LineWidth', 1.5);
lgd = legend([h_sec, h_top, h_nodes], {'Cylinder cross-sections', 'Top line of the cylinder', 'Nodal points'}, ...
    'Location', 'southwest');
set(lgd, 'Color', [0.2 0.2 0.2], 'TextColor', 'w');

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

% If no segment intersection found, check if endpoints are close (for base nodes)
endpoints1 = [points1(1,:); points1(end,:)];
endpoints2 = [points2(1,:); points2(end,:)];

for e1 = 1:2
    for e2 = 1:2
        dist = norm(endpoints1(e1,:) - endpoints2(e2,:));
        if dist < 50  % Tolerance in pixels for near-miss endpoints
            pt = (endpoints1(e1,:) + endpoints2(e2,:)) / 2;
            return;
        end
    end
end
end

function nodal_points = get_nodal_points(arcs_A, arcs_B, H)
% Find intersection between all pairs of arcs, optionally transforming them first
nodal_points = [];
for i = 1:length(arcs_A)
    pts_A = arcs_A{i};
    if nargin > 2 && ~isempty(H)
        pts_A_hom = [pts_A, ones(size(pts_A, 1), 1)]';
        pts_A_rect = H * pts_A_hom;
        pts_A = (pts_A_rect(1:2, :) ./ pts_A_rect(3, :))';
    end

    for j = 1:length(arcs_B)
        pts_B = arcs_B{j};
        if nargin > 2 && ~isempty(H)
            pts_B_hom = [pts_B, ones(size(pts_B, 1), 1)]';
            pts_B_rect = H * pts_B_hom;
            pts_B = (pts_B_rect(1:2, :) ./ pts_B_rect(3, :))';
        end

        pt = find_intersection(pts_A, pts_B);
        if ~isempty(pt)
            nodal_points = [nodal_points; pt, i, j];
        end
    end
end
end

function lambda = intersect_ray_cylinder(ray, P_axis, axis_dir, R)
% Find intersection of ray with cylinder
% Solves: ||cross(lambda*ray - P_axis, axis_dir)|| = R
% Returns NaN if the ray misses the cylinder.

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
    lambda = NaN;
    return;
end

% The camera is outside the cylinder and looks at the vault from below: the ray
% enters through the open lower half, so the vault is the farther intersection.
lambda = (-B + sqrt(discriminant)) / (2*A);
end

function line_hom = get_normalized_line(p1, p2, cx, cy, scale)
% Computes a homogeneous line representation in normalized coordinates.
p1_norm = [(p1(1)-cx)/scale, (p1(2)-cy)/scale, 1];
p2_norm = [(p2(1)-cx)/scale, (p2(2)-cy)/scale, 1];
line_hom = cross(p1_norm, p2_norm)';
end

function tf = crosses_line(points, l)
% True if the polyline has points on both sides of the homogeneous line l.
s = [points, ones(size(points, 1), 1)] * l(:);
tf = any(s > 0) && any(s < 0);
end

function pts = smooth_arc(points)
% Smooths the picked points of an arc with a degree-4 polynomial in chord length
% and resamples it densely.
s = [0; cumsum(sqrt(sum(diff(points).^2, 2)))];
s = s / s(end);
t = linspace(0, 1, 160)';
pts = [polyval(polyfit(s, points(:,1), 4), t), polyval(polyfit(s, points(:,2), 4), t)];
end

function mid = mirror_midpoints(arc_a, arc_b, vp)
% For every point p of arc_a, the line through p and the vanishing point vp is
% intersected with arc_b at q. Returns the image of the 3D midpoint of the two
% points, which is the harmonic conjugate of vp with respect to p and q.
mid = [];
for i = 1:size(arc_a, 1)
    p = arc_a(i, :);
    dir_p = vp - p;
    for j = 1:size(arc_b, 1) - 1
        b1 = arc_b(j, :);
        dir_b = arc_b(j+1, :) - b1;
        denom = dir_p(1)*dir_b(2) - dir_p(2)*dir_b(1);
        if abs(denom) < 1e-12, continue; end
        s = ((b1(1)-p(1))*dir_p(2) - (b1(2)-p(2))*dir_p(1)) / denom;
        if s < 0 || s > 1, continue; end
        q = b1 + s * dir_b;

        pq = q - p;
        t_vp = dot(vp - p, pq) / dot(pq, pq);
        mid = [mid; p + t_vp / (2*t_vp - 1) * pq];
    end
end
end

function lambda = fit_aspect_ratio(sections, H)
% Joint least-squares fit of x^2 + lambda^2 y^2 + a_k x + b_k y + c_k = 0 over all
% cross-sections, after mapping their points with H. Returns NaN if no ellipse fits.
n = numel(sections);
A = [];
for k = 1:n
    p = [sections(k).pts, ones(size(sections(k).pts, 1), 1)] * H';
    p = p(:, 1:2) ./ p(:, 3);
    rows_k = zeros(size(p, 1), 2 + 3*n);
    rows_k(:, 1:2) = p.^2;
    rows_k(:, 2 + 3*(k-1) + (1:3)) = [p, ones(size(p, 1), 1)];
    A = [A; rows_k];
end
col_scale = vecnorm(A);
[~, ~, V] = svd(A ./ col_scale, 'econ');
sol = V(:, end) ./ col_scale(:);
if sol(2) / sol(1) > 0
    lambda = sqrt(sol(2) / sol(1));
else
    lambda = NaN;
end
end

function row = iac_row(u, v)
% Coefficients of u' * omega * v for the zero-skew omega = [w1 0 w3; 0 w2 w4; w3 w4 w5].
row = [u(1)*v(1), u(2)*v(2), u(1)*v(3) + u(3)*v(1), u(2)*v(3) + u(3)*v(2), u(3)*v(3)];
end

function [center, radius] = fit_circle(x, y)
% Circle fit: algebraic start, then Gauss-Newton on the geometric distance.
abc = [x, y, ones(size(x))] \ -(x.^2 + y.^2);
p = [-abc(1)/2; -abc(2)/2; sqrt(abc(1)^2/4 + abc(2)^2/4 - abc(3))];
for iter = 1:20
    r = sqrt((x - p(1)).^2 + (y - p(2)).^2);
    J = [-(x - p(1)) ./ r, -(y - p(2)) ./ r, -ones(size(x))];
    p = p - J \ (r - p(3));
end
center = p(1:2)';
radius = p(3);
end

function [z, theta] = cylinder_coords(points, P_axis, axis_dir, vert_dir, trans_dir)
% Axial coordinate and angle around the axis (zero at the top) of 3D points.
rel = points - P_axis(:)';
z = rel * axis_dir(:);
theta = atan2(rel * trans_dir(:), rel * vert_dir(:));
end

function x = project_points(K, X, cols, rows)
% Projects 3D points (rows of X) to pixels. Points behind the camera or outside
% the image become NaN so that plotted curves break there.
x_hom = (K * X')';
x = x_hom(:, 1:2) ./ x_hom(:, 3);
outside = x_hom(:, 3) <= 0 | x(:,1) < 1 | x(:,1) > cols | x(:,2) < 1 | x(:,2) > rows;
x(outside, :) = NaN;
end

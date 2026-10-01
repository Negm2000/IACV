function [H_final, img_rect, out_size] = compute_metric_rectification(img, l_inf_perp, v_vert, v_trans)
% COMPUTE_METRIC_RECTIFICATION - Compute metric rectification homography
%   Returns the final homography, rectified image, and output size.

[rows, cols, ~] = size(img);

% Step 1: Affine rectification - map vanishing line to infinity
l_norm = l_inf_perp(:) / l_inf_perp(3);
H_aff = [1, 0, 0; 0, 1, 0; l_norm(1), l_norm(2), 1];

% Step 2: Metric rectification using known perpendicular directions
v_vert_aff = H_aff * v_vert(:);
v_trans_aff = H_aff * v_trans(:);

dir_vert = v_vert_aff(1:2) / norm(v_vert_aff(1:2));
dir_trans = v_trans_aff(1:2) / norm(v_trans_aff(1:2));

M = [dir_vert, dir_trans];
S_metric = [0, 1; 1, 0] * inv(M);

H_metric = [S_metric, [0; 0]; 0, 0, 1];
H_rect = H_metric * H_aff;

% Determine output view by sampling points
[xx, yy] = meshgrid(linspace(1, cols, 40), linspace(1, rows, 40));
pts_grid = [xx(:), yy(:), ones(numel(xx), 1)]';

grid_vals = H_aff(3,:) * pts_grid;
side = sign(mean(grid_vals));
if side == 0, side = 1; end

v_range = abs(max(grid_vals) - min(grid_vals));
margin = 0.15 * v_range;
mask = (sign(grid_vals) == side) & (abs(grid_vals) > margin);
if sum(mask) < 10, mask = (sign(grid_vals) == side); end

pts_trans = H_rect * pts_grid(:, mask);
pts_trans = pts_trans(1:2, :) ./ pts_trans(3, :);

% Bounding box
min_x = min(pts_trans(1,:)); max_x = max(pts_trans(1,:));
min_y = min(pts_trans(2,:)); max_y = max(pts_trans(2,:));
w_rect = max_x - min_x; h_rect_val = max_y - min_y;

if h_rect_val > 5 * w_rect, h_rect_val = 5 * w_rect; end
if w_rect > 5 * h_rect_val, w_rect = 5 * h_rect_val; end

% Scale and translate
target_w = 2000;
scale_factor = target_w / w_rect;
T = [scale_factor, 0, -scale_factor*min_x + 1; 0, scale_factor, -scale_factor*min_y + 1; 0, 0, 1];
H_final = T * H_rect;

out_size = [ceil(scale_factor * h_rect_val), target_w];
img_rect = imwarp(img, projective2d(H_final'), 'OutputView', imref2d(out_size));
end

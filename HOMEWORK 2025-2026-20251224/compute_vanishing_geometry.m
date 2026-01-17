function [v_vert, v_axis, v_trans, v_vert_n, v_axis_n, v_trans_n, l_inf_perp] = compute_vanishing_geometry(img, lines_v, lines_axis, line_apical, lines_trans)
% COMPUTE_VANISHING_GEOMETRY - Compute vanishing points and vanishing line
%   Returns both pixel-coordinate and normalized vanishing points.

[rows, cols, ~] = size(img);
cx = cols / 2;
cy = rows / 2;
scale = max(cx, cy);

T_n2p = [scale 0 cx; 0 scale cy; 0 0 1];

get_line_n = @(p1, p2) cross([(p1(1)-cx)/scale, (p1(2)-cy)/scale, 1], ...
    [(p2(1)-cx)/scale, (p2(2)-cy)/scale, 1])';

% Compute normalized vanishing points
v_vert_n = [];
if ~isempty(lines_v)
    lines_hom = {};
    for i=1:length(lines_v)
        lines_hom{i} = get_line_n(lines_v(i).p1, lines_v(i).p2);
    end
    v_vert_n = compute_vanishing_point(lines_hom);
end

v_axis_n = [];
lines_hom = {};
for i=1:length(lines_axis)
    lines_hom{end+1} = get_line_n(lines_axis(i).p1, lines_axis(i).p2);
end
if ~isempty(line_apical)
    lines_hom{end+1} = get_line_n(line_apical.p1, line_apical.p2);
end
if length(lines_hom) >= 2
    v_axis_n = compute_vanishing_point(lines_hom);
end

v_trans_n = [];
if ~isempty(lines_trans)
    lines_hom = {};
    for i=1:length(lines_trans)
        lines_hom{i} = get_line_n(lines_trans(i).p1, lines_trans(i).p2);
    end
    v_trans_n = compute_vanishing_point(lines_hom);
end

% Denormalize to pixel coordinates
v_vert = []; v_axis = []; v_trans = [];
if ~isempty(v_vert_n)
    v_vert = (T_n2p * v_vert_n(:))';
    v_vert = v_vert / v_vert(3);
end
if ~isempty(v_axis_n)
    v_axis = (T_n2p * v_axis_n(:))';
    v_axis = v_axis / v_axis(3);
end
if ~isempty(v_trans_n)
    v_trans = (T_n2p * v_trans_n(:))';
    v_trans = v_trans / v_trans(3);
end

% Compute vanishing line
l_inf_perp = [];
if ~isempty(v_vert) && ~isempty(v_trans)
    l_inf_perp = cross(v_vert, v_trans);
    l_inf_perp = l_inf_perp / norm(l_inf_perp(1:2));
end
end

function vp = compute_vanishing_point(lines_hom)
n = length(lines_hom);
if n < 2
    error('Need at least 2 lines.');
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

function nodal_points = find_nodal_points_geometric(arcs_A, arcs_B, H_R, v_axis, v_vert, line_apical)
% FIND_NODAL_POINTS_GEOMETRIC - Robust intersection finding with Index Alignment
%
% Usage:
%   nodal_points = find_nodal_points_geometric(arcs_A, arcs_B, H_R, v_axis, v_vert, line_apical)
%
% Features:
%   1. Finds actual geometric intersections of arc segments.
%   2. Uses 'line_apical' to detect if arc families are shifted (e.g., A1 matches B2).
%   3. Returns 'nodal_points' with LOGICAL indices, so N_ii is always Apical.

nodal_points = [];

if isempty(arcs_A) || isempty(arcs_B)
    fprintf('No arcs provided.\n');
    return;
end

if nargin < 3 || isempty(H_R), H_R = eye(3); end
% H_R is used here only if we wanted to filter by rectified geometry,
% but pure image-space intersection is often more robust for raw feature extraction.
H_R_inv = inv(H_R);

fprintf('Finding nodal points (Geometric Intersection)...\n');

%% 1. Transform arcs to Rectified Space for cleaner intersection checking
arcs_A_rect = cell(size(arcs_A));
arcs_B_rect = cell(size(arcs_B));

for i = 1:length(arcs_A)
    pts = arcs_A{i};
    pts_hom = [pts, ones(size(pts,1), 1)]';
    pts_rect = H_R * pts_hom;
    pts_rect = pts_rect(1:2, :) ./ pts_rect(3, :);
    arcs_A_rect{i} = pts_rect';
end

for i = 1:length(arcs_B)
    pts = arcs_B{i};
    pts_hom = [pts, ones(size(pts,1), 1)]';
    pts_rect = H_R * pts_hom;
    pts_rect = pts_rect(1:2, :) ./ pts_rect(3, :);
    arcs_B_rect{i} = pts_rect';
end

%% 2. Find ALL raw intersections
raw_crossings = []; % [x_img, y_img, index_A, index_B]

for i = 1:length(arcs_A_rect)
    pts_A = arcs_A_rect{i};
    for j = 1:length(arcs_B_rect)
        pts_B = arcs_B_rect{j};

        [pt_rect, found] = find_polyline_intersection(pts_A, pts_B);

        if found
            % Transform back to original image
            pt_orig = H_R_inv * [pt_rect(:); 1];
            pt_orig = pt_orig(1:2) / pt_orig(3);

            % Basic boundary check
            if pt_orig(1) > 0 && pt_orig(2) > 0 && pt_orig(1) < 10000
                raw_crossings = [raw_crossings; pt_orig', i, j];
            end
        end
    end
end

if isempty(raw_crossings)
    fprintf('No intersections found between arcs.\n');
    return;
end

fprintf('  Found %d raw intersections.\n', size(raw_crossings, 1));

%% 3. Auto-Align Indices using Apical Line
% Problem: User might have drawn A1, A2... but B2, B3... (skipping B1).
% Result: A1 intersects B2 (physically side node) but code thinks it's N(1,1) if shifted.
% Solution: The set of intersections that lies on 'line_apical' is the TRUE diagonal.

shift_detected = 0;

if nargin >= 6 && ~isempty(line_apical)
    fprintf('  Validating against Apical Line (Yellow)...\n');

    % Define the apical line function distance(pt)
    p1 = line_apical.p1; p2 = line_apical.p2;
    % Line eq: ax + by + c = 0
    a = p1(2) - p2(2);
    b = p2(1) - p1(1);
    c = -a*p1(1) - b*p1(2);
    norm_ab = sqrt(a^2 + b^2);
    dist_func = @(pt) abs(a*pt(1) + b*pt(2) + c) / norm_ab;

    % Check potential shifts k, where B_index = A_index + k
    % We search for the k that maximizes the number of points on the yellow line
    possible_shifts = -5:5;
    best_score = -1;
    best_k = 0;

    for k = possible_shifts
        % Get points where j - i = k
        candidates = raw_crossings(raw_crossings(:,4) - raw_crossings(:,3) == k, :);

        if isempty(candidates), continue; end

        % Count how many are close to apical line (e.g., within 20 pixels)
        dists = arrayfun(@(r) dist_func(candidates(r, 1:2)), 1:size(candidates,1));
        inliers = sum(dists < 30); % 30px tolerance

        if inliers > best_score
            best_score = inliers;
            best_k = k;
        end
    end

    if best_score > 0
        shift_detected = best_k;
        fprintf('  DETECTED INDEX SHIFT: k = %d (Best alignment with yellow line).\n', shift_detected);
        if shift_detected ~= 0
            fprintf('  -> Correcting: Physical N(i, i+%d) re-labeled as Apical N(i,i).\n', shift_detected);
        end
    else
        fprintf('  Warning: Could not align nodes with apical line. Assuming no shift.\n');
    end
else
    fprintf('  No apical line provided. Assuming indices are correct (k=0).\n');
end

%% 4. Build Final Output
% We adjust B-index so that N(i,i) is the apical one.
% Logical_B = Raw_B - shift
% Example: If shift=1 (A1 intersects B2 at apex), we want that point to be N(1,1).
% So we output B_index = 2 - 1 = 1.

for k = 1:size(raw_crossings, 1)
    pt = raw_crossings(k, 1:2);
    idx_A = raw_crossings(k, 3);
    idx_B_raw = raw_crossings(k, 4);

    idx_B_logical = idx_B_raw - shift_detected;

    % Only keep if the logical index is valid?
    % Actually, keep all, even if indices go negative (user debug info),
    % but strictly speaking 'reconstruct_3d' needs valid matches.

    nodal_points = [nodal_points; pt, idx_A, idx_B_logical];

    type = 'Side';
    if idx_A == idx_B_logical, type = 'APICAL'; end

    % Debug print for Apical ones
    if strcmp(type, 'APICAL')
        fprintf('  %s Node N(%d,%d) at (%.0f,%.0f)\n', type, idx_A, idx_B_logical, pt(1), pt(2));
    end
end

end

%% ========== HELPER FUNCTION ==========

function [pt, found] = find_polyline_intersection(pts1, pts2)
% Standard segment-to-segment intersection
pt = []; found = false;
n1 = size(pts1, 1); n2 = size(pts2, 1);
if n1 < 2 || n2 < 2, return; end

for i = 1:n1-1
    p1a = pts1(i, :); p1b = pts1(i+1, :);
    for j = 1:n2-1
        p2a = pts2(j, :); p2b = pts2(j+1, :);
        d1 = p1b - p1a; d2 = p2b - p2a; d12 = p1a - p2a;
        denom = d1(1)*d2(2) - d1(2)*d2(1);
        if abs(denom) < 1e-10, continue; end
        t1 = (d2(1)*d12(2) - d2(2)*d12(1)) / denom;
        t2 = (d1(1)*d12(2) - d1(2)*d12(1)) / denom;
        if t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1
            pt = p1a + t1 * d1; found = true; return;
        end
    end
end
end

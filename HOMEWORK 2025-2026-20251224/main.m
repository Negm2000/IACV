% IACV Homework 2025-2026
% 3D Reconstruction of a Cylindric Vault
% Author: [Your Name]
% Date: 2026-01-03

clear; close all; clc;

%% Setup and Parameters
imageFile = 'San Maurizio.jpg';
featureFile = 'features_fixed.mat';

%% 1. Load Image
if exist(imageFile, 'file')
    img = imread(imageFile);
    % figure(1); imshow(img); title('Input Image: San Maurizio');
    % hold on;
else
    error('Image file %s not found. Please make sure it is in the current directory.', imageFile);
end

[rows, cols, ~] = size(img);

%% 2. Feature Extraction (Smart Loading Version)
% This section checks what you HAVE vs what you NEED.
% Initialize empty/default values
if ~exist('lines_v', 'var'), lines_v = []; end
if ~exist('lines_axis', 'var'), lines_axis = []; end
if ~exist('line_apical', 'var'), line_apical = []; end
if ~exist('lines_trans', 'var'), lines_trans = []; end
if ~exist('arcs_A', 'var'), arcs_A = {}; end
if ~exist('arcs_B', 'var'), arcs_B = {}; end

needs_save = false;

% 1. Try to load existing data
if exist(featureFile, 'file')
    fprintf('Loading existing features from %s...\n', featureFile);
    load(featureFile);
else
    fprintf('No feature file found. Starting fresh.\n');
end

% 2. Fast Re-pick Selection
repickList = [];
if exist(featureFile, 'file')
    prompt = {['Select features to RE-PICK (type numbers e.g. "1 3 5" or leave empty to KEEP ALL):', newline, ...
        '1. Vertical (Black)', newline, ...
        '2. Axis-Parallel (Green)', newline, ...
        '3. Apical (Yellow)', newline, ...
        '4. Transversal (White)', newline, ...
        '5. Arcs Family A (Cyan)', newline, ...
        '6. Arcs Family B (Magenta)']};
    answer = inputdlg(prompt, 'Refine Features', [1 100]);

    if ~isempty(answer) && ~isempty(answer{1})
        repickList = str2num(answer{1}); % Parse string like "1 3 5" to [1 3 5]
    end
end

% --- Vertical Lines ---
if isempty(lines_v) || ismember(1, repickList)
    fprintf('1. Select VERTICAL lines (Black).\n');
    lines_v = select_lines(img, 'Select Vertical Lines (Black). Enter to finish.', 'k');
    needs_save = true;
end

% --- Axis Parallel Lines ---
if isempty(lines_axis) || ismember(2, repickList)
    fprintf('2. Select AXIS-PARALLEL lines (Green).\n');
    lines_axis = select_lines(img, 'Select Axis-Parallel Lines (Green). Enter to finish.', 'g');
    needs_save = true;
end

% --- Apical Line ---
if isempty(line_apical) || ismember(3, repickList)
    fprintf('3. Select THE APICAL line (Yellow).\n');
    temp_apical = select_lines(img, 'Select THE Apical Line (Yellow). Select ONE and Enter.', 'y');
    if ~isempty(temp_apical), line_apical = temp_apical(1); end
    needs_save = true;
end

% --- Transversal Lines ---
if isempty(lines_trans) || ismember(4, repickList)
    fprintf('4. Select TRANSVERSAL lines (White).\n');
    lines_trans = select_lines(img, 'Select Transversal Lines (White). Enter to finish.', 'w');
    needs_save = true;
end

% --- Arcs Family A ---
if isempty(arcs_A) || ismember(5, repickList)
    fprintf('5. Select points for DIAGONAL ARCS (Family A - Cyan).\n');
    arcs_A = select_arcs(img, 'Select Family A Arcs (Up-Right /). Enter on empty to Finish');
    needs_save = true;
end

% --- Arcs Family B ---
if isempty(arcs_B) || ismember(6, repickList)
    fprintf('6. Select points for DIAGONAL ARCS (Family B - Magenta).\n');
    arcs_B = select_arcs(img, 'Select Family B Arcs (Up-Left \). Enter twice to Finish.');
    needs_save = true;
end

% 3. Save if we added anything new
if needs_save
    save(featureFile, 'lines_v', 'lines_axis', 'line_apical', 'lines_trans', 'arcs_A', 'arcs_B');
    fprintf('Features updated and saved to %s\n', featureFile);
end

% Visualizer of everything
figure(1); imshow(img); title('All Features Loaded'); hold on;
if ~isempty(lines_v), for i=1:length(lines_v), plot([lines_v(i).p1(1) lines_v(i).p2(1)], [lines_v(i).p1(2) lines_v(i).p2(2)], 'k-', 'LineWidth', 2); end; end
if ~isempty(lines_axis), for i=1:length(lines_axis), plot([lines_axis(i).p1(1) lines_axis(i).p2(1)], [lines_axis(i).p1(2) lines_axis(i).p2(2)], 'g-', 'LineWidth', 2); end; end
if ~isempty(line_apical), plot([line_apical.p1(1) line_apical.p2(1)], [line_apical.p1(2) line_apical.p2(2)], 'y-', 'LineWidth', 3); end
if ~isempty(lines_trans), for i=1:length(lines_trans), plot([lines_trans(i).p1(1) lines_trans(i).p2(1)], [lines_trans(i).p1(2) lines_trans(i).p2(2)], 'w-', 'LineWidth', 2); end; end
for i=1:length(arcs_A), pts=arcs_A{i}; plot(pts(:,1), pts(:,2), 'c.-', 'MarkerSize', 10); end
for i=1:length(arcs_B), pts=arcs_B{i}; plot(pts(:,1), pts(:,2), 'm.-', 'MarkerSize', 10); end

%% 3. Vanishing Points and Vanishing Line
fprintf('Computing Vanishing Points...\n');

% Helper function to get homogeneous line from two points
get_line = @(p1, p2) cross([p1, 1], [p2, 1])';

% --- v_vert (Black) ---
if ~isempty(lines_v)
    lines_hom = {};
    for i=1:length(lines_v), lines_hom{i} = get_line(lines_v(i).p1, lines_v(i).p2); end
    v_vert = compute_vanishing_point(lines_hom);
    fprintf('v_vert  = [%.4f, %.4f, %.4f]\n', v_vert);
else
    v_vert = [];
end

% --- v_axis (Green + Yellow Combined) ---
lines_hom = {};
for i=1:length(lines_axis), lines_hom{end+1} = get_line(lines_axis(i).p1, lines_axis(i).p2); end
if ~isempty(line_apical), lines_hom{end+1} = get_line(line_apical.p1, line_apical.p2); end

if length(lines_hom) >= 2
    v_axis = compute_vanishing_point(lines_hom);
    fprintf('v_axis  = [%.4f, %.4f, %.4f]\n', v_axis);
else
    v_axis = [];
end

% --- v_trans (White) ---
if ~isempty(lines_trans)
    lines_hom = {};
    for i=1:length(lines_trans), lines_hom{i} = get_line(lines_trans(i).p1, lines_trans(i).p2); end
    v_trans = compute_vanishing_point(lines_hom);
    fprintf('v_trans = [%.4f, %.4f, %.4f]\n', v_trans);
else
    v_trans = [];
end

% Compute Vanishing Line of planes perpendicular to the cylinder axis
if ~isempty(v_vert) && ~isempty(v_trans)
    l_inf_perp = cross(v_vert, v_trans);
    l_inf_perp = l_inf_perp / norm(l_inf_perp(1:2));
    fprintf('Vanishing line l_inf_perp = [%.4f, %.4f, %.4f]\n', l_inf_perp);
else
    l_inf_perp = [];
end

% Visualize vanishing points and line extensions
if ~isempty(v_vert) || ~isempty(v_axis) || ~isempty(v_trans)
    figure(2); clf; imshow(img); title('Vanishing Points & Line Extensions'); hold on;

    % Plotting limits calculation
    all_pts = [1, 1; cols, rows]; % Start with image corners

    % --- v_vert (Black) ---
    if ~isempty(v_vert) && abs(v_vert(3)) > 1e-9
        vp = v_vert(1:2)/v_vert(3); all_pts = [all_pts; vp];
        plot(vp(1), vp(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
        text(vp(1), vp(2), '  v_{vert}', 'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
        for i=1:length(lines_v)
            mid = (lines_v(i).p1 + lines_v(i).p2) / 2;
            plot([mid(1), vp(1)], [mid(2), vp(2)], 'k:', 'LineWidth', 0.5);
            % Intersect with l_inf_perp (v_vert is already on it, but for consistency)
            if ~isempty(l_inf_perp)
                L = cross([lines_v(i).p1 1], [lines_v(i).p2 1]);
                pt_inf = cross(L, l_inf_perp);
                if abs(pt_inf(3)) > 1e-9
                    pt_inf = pt_inf(1:2)/pt_inf(3); all_pts = [all_pts; pt_inf];
                    plot([mid(1), pt_inf(1)], [mid(2), pt_inf(2)], 'k:', 'LineWidth', 0.5);
                end
            end
        end
    end

    % --- v_axis (Green) ---
    if ~isempty(v_axis) && abs(v_axis(3)) > 1e-9
        vp = v_axis(1:2)/v_axis(3); all_pts = [all_pts; vp];
        plot(vp(1), vp(2), 'go', 'MarkerSize', 12, 'LineWidth', 2);
        text(vp(1), vp(2), '  v_{axis}', 'Color', 'g', 'FontSize', 10, 'FontWeight', 'bold');
        for i=1:length(lines_axis)
            mid = (lines_axis(i).p1 + lines_axis(i).p2) / 2;
            plot([mid(1), vp(1)], [mid(2), vp(2)], 'g:', 'LineWidth', 0.5);
            % Extend to l_inf_perp
            if ~isempty(l_inf_perp)
                L = cross([lines_axis(i).p1 1], [lines_axis(i).p2 1]);
                pt_inf = cross(L, l_inf_perp);
                if abs(pt_inf(3)) > 1e-9
                    pt_inf = pt_inf(1:2)/pt_inf(3); all_pts = [all_pts; pt_inf];
                    plot([mid(1), pt_inf(1)], [mid(2), pt_inf(2)], 'g:', 'LineWidth', 0.5);
                end
            end
        end
        if ~isempty(line_apical)
            mid = (line_apical.p1 + line_apical.p2) / 2;
            plot([mid(1), vp(1)], [mid(2), vp(2)], 'y:', 'LineWidth', 0.8);
            if ~isempty(l_inf_perp)
                L = cross([line_apical.p1 1], [line_apical.p2 1]);
                pt_inf = cross(L, l_inf_perp);
                if abs(pt_inf(3)) > 1e-9
                    pt_inf = pt_inf(1:2)/pt_inf(3); all_pts = [all_pts; pt_inf];
                    plot([mid(1), pt_inf(1)], [mid(2), pt_inf(2)], 'y:', 'LineWidth', 0.8);
                end
            end
        end
    end

    % --- v_trans (White) ---
    if ~isempty(v_trans) && abs(v_trans(3)) > 1e-9
        vp = v_trans(1:2)/v_trans(3); all_pts = [all_pts; vp];
        plot(vp(1), vp(2), 'ro', 'MarkerSize', 12, 'LineWidth', 2); % Red for contrast on white
        text(vp(1), vp(2), '  v_{trans}', 'Color', [0.8 0 0], 'FontSize', 10, 'FontWeight', 'bold');
        for i=1:length(lines_trans)
            mid = (lines_trans(i).p1 + lines_trans(i).p2) / 2;
            plot([mid(1), vp(1)], [mid(2), vp(2)], 'w:', 'LineWidth', 0.5);
            % Extend to l_inf_perp
            if ~isempty(l_inf_perp)
                L = cross([lines_trans(i).p1 1], [lines_trans(i).p2 1]);
                pt_inf = cross(L, l_inf_perp);
                if abs(pt_inf(3)) > 1e-9
                    pt_inf = pt_inf(1:2)/pt_inf(3); all_pts = [all_pts; pt_inf];
                    plot([mid(1), pt_inf(1)], [mid(2), pt_inf(2)], 'w:', 'LineWidth', 0.5);
                end
            end
        end
    end

    % Adjust Axis to show image + points with slight padding
    x_min = min(all_pts(:,1)); x_max = max(all_pts(:,1));
    y_min = min(all_pts(:,2)); y_max = max(all_pts(:,2));
    pad_x = abs(x_max - x_min) * 0.1 + 100; % Larger padding
    pad_y = abs(y_max - y_min) * 0.1 + 100;
    curr_axis = [x_min-pad_x, x_max+pad_x, y_min-pad_y, y_max+pad_y];

    % --- FIX: Ensure axis is visible and grid is on top ---
    set(gca, 'Visible', 'on'); % Show coordinates and axis lines
    axis(curr_axis);
    set(gca, 'Layer', 'top');  % Draw grid OVER features if they overlap
    grid on; grid minor;
    set(gca, 'Color', [0.9 0.9 0.9]); % Grey background for the "infinite" space

    % --- Vanishing Line l_inf_perp ---
    if ~isempty(l_inf_perp)
        % Draw the line across the current axis
        % Line eq: ax + by + c = 0
        a = l_inf_perp(1); b = l_inf_perp(2); c = l_inf_perp(3);

        % Sample two points at the edge of the current x-limits
        xl = [curr_axis(1), curr_axis(2)];
        if abs(b) > 1e-9
            yl = -(a*xl + c)/b;
            plot(xl, yl, 'b-', 'LineWidth', 2.5);
            text(xl(2), yl(2), '  l_{\infty\perp}', 'Color', 'b', 'FontSize', 14, 'FontWeight', 'bold');
        else
            % Vertical line case
            x_vert = -c/a;
            plot([x_vert, x_vert], [curr_axis(3), curr_axis(4)], 'b-', 'LineWidth', 2.5);
            text(x_vert, curr_axis(4), '  l_{\infty\perp}', 'Color', 'b', 'FontSize', 14, 'FontWeight', 'bold');
        end
    end
end

%% 5. Rectification (Metric Rectification of Transversal Plane)
% Goal: Rectify the plane perpendicular to the cylinder axis.
% This plane contains v_vert and v_trans, with vanishing line l_inf_perp.
fprintf('Computing Rectification...\n');

H_R = [];
img_rect = [];

if ~isempty(l_inf_perp) && ~isempty(v_vert) && ~isempty(v_trans)
    % --- Step 1: Affine Rectification (map l_inf_perp to [0,0,1]) ---
    % H_aff maps l_inf_perp (the vanishing line) to the ideal line at infinity.
    % Theory: If l' = H^(-T) * l, then for l'=[0;0;1], l = H^T * [0;0;1] = 3rd col of H^T = 3rd row of H.
    % So we need H such that H(3,:) is proportional to l_inf_perp.

    l = l_inf_perp(:); % Ensure column
    l = l / norm(l(1:2)); % Normalize for numerical stability

    % Construct an affine homography H_aff
    % H_aff = [1 0 0; 0 1 0; l1 l2 l3] where l_inf_perp = [l1 l2 l3]'
    H_aff = [1 0 0; 0 1 0; l(1) l(2) l(3)];

    % --- Step 2: Metric Rectification (restore angles) ---
    % After H_aff, the vanishing points v_vert and v_trans become ideal points (at infinity).
    % In the affine-rectified image, they define directions.
    % We use their orthogonality to find the metric correction.

    % Transform vanishing points to affine-rectified space
    v_vert_aff = H_aff * v_vert(:);
    v_trans_aff = H_aff * v_trans(:);

    fprintf('v_vert_aff = [%.4f, %.4f, %.4f]\n', v_vert_aff);
    fprintf('v_trans_aff = [%.4f, %.4f, %.4f]\n', v_trans_aff);

    % These should be ideal points (w ≈ 0). Check and warn if not.
    if abs(v_vert_aff(3)) > 1e-3 || abs(v_trans_aff(3)) > 1e-3
        warning('Affine transformation did not map VPs to infinity as expected.');
    end

    % Extract directions (x, y components only)
    d1 = v_vert_aff(1:2); d1 = d1 / norm(d1);
    d2 = v_trans_aff(1:2); d2 = d2 / norm(d2);

    fprintf('d1 = [%.4f, %.4f], d2 = [%.4f, %.4f]\n', d1(1), d1(2), d2(1), d2(2));
    fprintf('Dot product (should be ~0 for orthogonal dirs): %.6f\n', dot(d1, d2));

    % Orthogonality constraint: d1' * S * d2 = 0, where S = [s1 s2; s2 s3].
    % We parameterize S symmetrically: s = [s1; s2; s3].
    % Constraint: s1*d1x*d2x + s2*(d1x*d2y + d1y*d2x) + s3*d1y*d2y = 0
    %
    % We have only ONE orthogonality constraint, but S has 3 unknowns (2 DoF up to scale).
    % Standard fix: Use the Image of Circular Points constraint.
    % For a simple metric rectification, we can assume the affine distortion is pure scaling:
    % S = [s 0; 0 1]. This gives: s*d1x*d2x + d1y*d2y = 0 => s = -d1y*d2y / (d1x*d2x).
    %
    % If s < 0, it means d1 and d2 are NOT orthogonal in the original scene (or VP error).
    % In this case, we use a rotation matrix to make them orthogonal.

    denom = d1(1)*d2(1);
    if abs(denom) > 1e-9
        s_ratio = -(d1(2)*d2(2)) / denom;
    else
        s_ratio = 1;
        warning('Metric rectification: denominator near zero, using s=1.');
    end

    fprintf('Computed s_ratio = %.6f\n', s_ratio);

    if s_ratio > 0
        % Standard case: apply scaling
        A_metric = [sqrt(s_ratio) 0; 0 1];
    else
        % s < 0 means the VPs are not orthogonal in the affine space.
        % Apply a rotation to align d1 with the y-axis, then scale.
        % This effectively rotates the image so that v_vert is truly vertical.

        theta = atan2(d1(1), d1(2)); % Angle of d1 from vertical (y-axis)
        R = [cos(theta) sin(theta); -sin(theta) cos(theta)]; % Rotate d1 to [0; 1]

        % After rotation, recompute d2 and the scaling
        d2_rot = R * d2;

        % Now d1 is [0; 1], so the orthogonality constraint simplifies:
        % For d1' * S * d2 = 0 with d1 = [0; 1]: s3 * d2y = 0.
        % If d2y ≠ 0, we have s3 = 0, which is degenerate.
        % Better approach: just use the rotation to correct the metric.

        fprintf('Using rotation-based metric correction (theta = %.2f deg).\n', rad2deg(theta));
        A_metric = R; % Use rotation as the metric correction
    end

    % The metric homography in the affine-rectified frame is:
    H_metric = [A_metric, [0;0]; 0 0 1];

    % --- Step 3: Combine and Apply ---
    H_R = H_metric * H_aff;

    fprintf('Rectification Homography H_R computed.\n');
    disp(H_R);

    % Apply to image
    tform = projective2d(H_R');
    img_rect = imwarp(img, tform, 'OutputView', imref2d(size(img)));

    % Visualize
    figure(3); clf;
    subplot(1,2,1); imshow(img); title('Original Image');
    subplot(1,2,2); imshow(img_rect); title('Rectified Image (Transversal Plane)');
else
    warning('Could not compute rectification: Missing l_inf_perp, v_vert, or v_trans.');
end

%% 6. Camera Calibration (K)
% Compute matrix K using orthogonality of {v_vert, v_axis, v_trans}.
% Theory: v_i' * omega * v_j = 0.
% K depends on four parameters: fx, fy, Uo, Vo (Zero skew).
% Thus omega = [w1 0 w2; 0 w3 w4; w2 w4 w5] (5 non-zero entries).
% Note: w1 != w3 because fx != fy.

fprintf('Computing Calibration Matrix K...\n');
% TODO: Solve linear system A*x = 0 for omega entries.
% K = inv(chol(omega));

%% 6. 3D Reconstruction
% 6.1 Back-projection: Ray direction d = K \ [u; v; 1]
% 6.2 Intersection with Cylinder Model: X^2 + Y^2 = R^2
% 6.3 Localize cylinder axis using apical points information.

fprintf('Performing 3D Reconstruction...\n');
% TODO: Implementation of ray-cylinder intersection logic.

%% Part 2 - Visualization

%% 7. Plotting
% Plot one of the curved arcs and show different views in 3D.

figure('Name', '3D Reconstruction');
% plot3(...)
grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');

disp('Processing Complete.');

%% Helper Functions for Feature Extraction

function lines = select_lines(img, promptStr, colorCode)
if nargin < 3, colorCode = 'y'; end % Default yellow
hFig = figure('Name', promptStr); imshow(img); hold on;
title({promptStr, 'MOUSE: Zoom/Pan with Toolbar tools', ...
    'ACTION: Toggle tools OFF to click/pick points', ...
    '[u]: Undo, [Enter]: Finish category'});

lines = struct('p1', {}, 'p2', {});
plotHandles = {}; % Track plot handles for undo
tempHandle = [];

while true
    % waitforbuttonpress allows toolbar tools to function!
    k = waitforbuttonpress;

    if k == 0 % Mouse click
        % CRITICAL: Only pick if NOT in Zoom or Pan mode
        zoomObj = zoom(hFig);
        panObj = pan(hFig);
        if strcmp(zoomObj.Enable, 'on') || strcmp(panObj.Enable, 'on')
            continue; % Let the tool handle it
        end

        currPt = get(gca, 'CurrentPoint');
        x = currPt(1,1);
        y = currPt(1,2);

        if isempty(tempHandle)
            % First point: Use high-precision '+' marker
            tempHandle = plot(x, y, [colorCode '+'], 'MarkerSize', 12, 'LineWidth', 1.5);
            p1 = [x, y];
            fprintf('Point 1 set at (%.2f, %.2f)\n', x, y);
        else
            % Second point
            p2 = [x, y];
            % --- FIX: Store both points in a SINGLE struct entry ---
            idx = length(lines) + 1;
            lines(idx).p1 = p1;
            lines(idx).p2 = p2;
            % --------------------------------------------------------
            delete(tempHandle);
            tempHandle = [];
            % Plot the line and endpoints
            hLine = plot([p1(1), p2(1)], [p1(2), p2(2)], [colorCode '-'], 'LineWidth', 2);
            hP1 = plot(p1(1), p1(2), [colorCode '+'], 'MarkerSize', 8);
            hP2 = plot(p2(1), p2(2), [colorCode '+'], 'MarkerSize', 8);
            plotHandles{end+1} = [hLine, hP1, hP2];
            fprintf('Line added.\n');
        end

    elseif k == 1 % Key press
        key = get(hFig, 'CurrentCharacter');
        if isempty(key), continue; end

        if key == 13 % Enter
            if isempty(tempHandle)
                break; % Finished
            else
                fprintf('Current line incomplete. Pick second point or [u] to cancel.\n');
            end
        elseif key == 'u' % Undo
            if ~isempty(tempHandle)
                delete(tempHandle);
                tempHandle = [];
                fprintf('Point 1 cancelled.\n');
            elseif ~isempty(lines)
                lines(end) = [];
                delete(plotHandles{end});
                plotHandles(end) = [];
                fprintf('Last line undone.\n');
            else
                fprintf('Nothing to undo.\n');
            end
        end
    end
end
close(hFig);
end

function arcs = select_arcs(img, promptStr)
% Returns cell array of point sets [x, y]
hFig = figure('Name', promptStr); imshow(img); hold on;
title({promptStr, 'MOUSE: Click to pick (ensure Zoom tool is OFF)', ...
    '[Backspace]: Undo point, [u]: Undo Arc, [Enter]: Save Arc/Finish'});

arcs = {};
arcPlots = {};
currentArcPts = [];
currentArcPlot = [];

while true
    title({promptStr, sprintf('Arc %d: %d points picked', length(arcs)+1, size(currentArcPts, 1))});

    k = waitforbuttonpress;
    if k == 0 % Click
        zoomObj = zoom(hFig); panObj = pan(hFig);
        if strcmp(zoomObj.Enable, 'on') || strcmp(panObj.Enable, 'on'), continue; end

        currPt = get(gca, 'CurrentPoint');
        x = currPt(1,1); y = currPt(1,2);
        currentArcPts = [currentArcPts; x, y];
        if ~isempty(currentArcPlot), delete(currentArcPlot); end
        currentArcPlot = plot(currentArcPts(:,1), currentArcPts(:,2), 'r.-', 'MarkerSize', 12, 'LineWidth', 1.5);
        fprintf('Point picked at (%.1f, %.1f)\n', x, y);

    elseif k == 1 % Key press
        key = get(hFig, 'CurrentCharacter');
        if isempty(key), continue; end

        if key == 13 % Enter
            if isempty(currentArcPts)
                break; % Exit
            else
                arcs{end+1} = currentArcPts;
                arcPlots{end+1} = currentArcPlot;
                currentArcPts = [];
                currentArcPlot = [];
                fprintf('Arc saved.\n');
            end
        elseif key == 'u' || key == 8 % 'u' or Backspace
            if ~isempty(currentArcPts)
                currentArcPts(end,:) = [];
                delete(currentArcPlot);
                if ~isempty(currentArcPts)
                    currentArcPlot = plot(currentArcPts(:,1), currentArcPts(:,2), 'r.-', 'MarkerSize', 12);
                else
                    currentArcPlot = [];
                end
                fprintf('Point undone.\n');
            elseif ~isempty(arcs)
                delete(arcPlots{end});
                arcs(end) = [];
                arcPlots(end) = [];
                fprintf('Whole arc removed.\n');
            end
        end
    end
end
close(hFig);
end

function vp = compute_vanishing_point(lines_hom)
% Compute vanishing point from a cell array of homogeneous lines.
% Uses least-squares intersection (SVD) for robustness.
% Each line is a 3x1 vector [a; b; c] representing ax + by + c = 0.

n = length(lines_hom);
if n < 2
    error('Need at least 2 lines to compute a vanishing point.');
end

% Build matrix A where each row is a line
A = zeros(n, 3);
for k = 1:n
    A(k, :) = lines_hom{k}';
end

% Solve A * vp = 0 using SVD (vp is the null space of A)
[~, ~, V] = svd(A);
vp = V(:, end); % Last column of V is the null space

% Normalize so that vp(3) = 1 if possible (for finite points)
if abs(vp(3)) > 1e-9
    vp = vp / vp(3);
end
vp = vp';
end

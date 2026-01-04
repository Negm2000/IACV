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
    figure(1); imshow(img); title('Input Image: San Maurizio');
    hold on;
else
    error('Image file %s not found. Please make sure it is in the current directory.', imageFile);
end

[rows, cols, ~] = size(img);

%% 2. Feature Extraction (Lines and Points)
% We need to extract:
% - Vertical lines (Black): giving v_vert
% - Axis-parallel lines (Green): generic parallel lines, giving v_axis
% - Apical line (Yellow): specifically connects highest nodal points N_ii
% - Transversal lines (White): connecting symmetric nodal points, giving v_trans
% - Points on diagonal arcs (for 3D reconstruction)

if exist(featureFile, 'file')
    fprintf('Loading features from %s...\n', featureFile);
    load(featureFile);

    % Visualize loaded features to verify
    figure(1); imshow(img); title('Loaded Features'); hold on;
    % Plot Vertical Lines (Black)
    if exist('lines_v', 'var')
        for i = 1:length(lines_v)
            plot([lines_v(i).p1(1), lines_v(i).p2(1)], [lines_v(i).p1(2), lines_v(i).p2(2)], 'k-', 'LineWidth', 2);
        end
    end
    % Plot Axis-Parallel Lines (Green)
    if exist('lines_axis', 'var')
        for i = 1:length(lines_axis)
            plot([lines_axis(i).p1(1), lines_axis(i).p2(1)], [lines_axis(i).p1(2), lines_axis(i).p2(2)], 'g-', 'LineWidth', 2);
        end
    end
    % Plot Apical Line (Yellow)
    if exist('line_apical', 'var') && ~isempty(line_apical)
        plot([line_apical.p1(1), line_apical.p2(1)], [line_apical.p1(2), line_apical.p2(2)], 'y-', 'LineWidth', 3);
    end
    % Plot Transversal Lines (White)
    if exist('lines_trans', 'var')
        for i = 1:length(lines_trans)
            plot([lines_trans(i).p1(1), lines_trans(i).p2(1)], [lines_trans(i).p1(2), lines_trans(i).p2(2)], 'w-', 'LineWidth', 2);
        end
    end
    % Plot Arcs (Points)
    if exist('arcs_A', 'var')
        colors = {'r.', 'm.', 'c.'};
        for i = 1:length(arcs_A)
            pts = arcs_A{i};
            plot(pts(:,1), pts(:,2), colors{mod(i-1, length(colors))+1}, 'MarkerSize', 10);
        end
    end

else
    fprintf('Feature file not found. Starting manual extraction...\n');

    fprintf('1. Select VERTICAL lines (Black in instructions).\n');
    lines_v = select_lines(img, 'Select Vertical Lines (Black). Press Enter when done.', 'k');

    fprintf('2. Select AXIS-PARALLEL lines (Green in instructions).\n');
    lines_axis = select_lines(img, 'Select Axis-Parallel Lines (Green). Press Enter when done.', 'g');

    fprintf('3. Select THE APICAL line (Yellow in instructions).\n');
    temp_apical = select_lines(img, 'Select THE Apical Line (Yellow). Select ONE and press Enter.', 'y');
    if ~isempty(temp_apical)
        line_apical = temp_apical(1);
    else
        line_apical = [];
    end

    fprintf('4. Select TRANSVERSAL lines (White in instructions).\n');
    lines_trans = select_lines(img, 'Select Transversal Lines (White). Press Enter when done.', 'w');

    fprintf('5. Select points for DIAGONAL ARCS.\n');
    arcs_A = select_arcs(img, 'Select Points for Diagonal Arcs. Press Enter after each arc.');

    save(featureFile, 'lines_v', 'lines_axis', 'line_apical', 'lines_trans', 'arcs_A');
    fprintf('Features saved to %s\n', featureFile);
end

%% 3. Vanishing Points and Vanishing Line
% 3.1 Compute Orthogonal Vanishing Points: v_vert, v_axis, v_trans
% 3.2 Find vanishing line l_inf of planes perpendicular to cylinder axis
%     (Planes perpendicular to axis are spanned by v_vert and v_trans)

fprintf('Computing Vanishing Points...\n');

% Helper function to get homogeneous line from two points
get_line = @(p1, p2) cross([p1, 1], [p2, 1])';

% Compute v_vert from vertical lines
if exist('lines_v', 'var') && length(lines_v) >= 2
    lines_v_hom = cell(1, length(lines_v));
    for idx = 1:length(lines_v)
        lines_v_hom{idx} = get_line(lines_v(idx).p1, lines_v(idx).p2);
    end
    v_vert = compute_vanishing_point(lines_v_hom);
    fprintf('v_vert = [%.4f, %.4f, %.4f]\n', v_vert);
else
    warning('Not enough vertical lines to compute v_vert.');
    v_vert = [];
end

% Compute v_axis from axis-parallel lines (Green)
if exist('lines_axis', 'var') && length(lines_axis) >= 2
    lines_axis_hom = cell(1, length(lines_axis));
    for idx = 1:length(lines_axis)
        lines_axis_hom{idx} = get_line(lines_axis(idx).p1, lines_axis(idx).p2);
    end
    v_axis = compute_vanishing_point(lines_axis_hom);
    fprintf('v_axis = [%.4f, %.4f, %.4f]\n', v_axis);
else
    warning('Not enough axis-parallel lines to compute v_axis.');
    v_axis = [];
end

% Compute v_trans from transversal lines (White)
if exist('lines_trans', 'var') && length(lines_trans) >= 2
    lines_trans_hom = cell(1, length(lines_trans));
    for idx = 1:length(lines_trans)
        lines_trans_hom{idx} = get_line(lines_trans(idx).p1, lines_trans(idx).p2);
    end
    v_trans = compute_vanishing_point(lines_trans_hom);
    fprintf('v_trans = [%.4f, %.4f, %.4f]\n', v_trans);
else
    warning('Not enough transversal lines to compute v_trans.');
    v_trans = [];
end

% Compute Vanishing Line of planes perpendicular to the cylinder axis
% These planes are spanned by Vertical and Transversal directions
if ~isempty(v_vert) && ~isempty(v_trans)
    l_inf_perp = cross(v_vert, v_trans);
    l_inf_perp = l_inf_perp / norm(l_inf_perp(1:2)); % Normalize
    fprintf('Vanishing line l_inf_perp = [%.4f, %.4f, %.4f]\n', l_inf_perp);
else
    warning('Cannot compute l_inf_perp without v_vert and v_trans.');
    l_inf_perp = [];
end

% Visualize vanishing points on the image
figure(2); imshow(img); title('Vanishing Points'); hold on;
if ~isempty(v_vert) && abs(v_vert(3)) > 1e-6
    plot(v_vert(1)/v_vert(3), v_vert(2)/v_vert(3), 'ko', 'MarkerSize', 15, 'LineWidth', 2);
    text(v_vert(1)/v_vert(3), v_vert(2)/v_vert(3), ' v_{vert}', 'Color', 'k', 'FontSize', 12);
end
if ~isempty(v_axis) && abs(v_axis(3)) > 1e-6
    plot(v_axis(1)/v_axis(3), v_axis(2)/v_axis(3), 'go', 'MarkerSize', 15, 'LineWidth', 2);
    text(v_axis(1)/v_axis(3), v_axis(2)/v_axis(3), ' v_{axis}', 'Color', 'g', 'FontSize', 12);
end
if ~isempty(v_trans) && abs(v_trans(3)) > 1e-6
    plot(v_trans(1)/v_trans(3), v_trans(2)/v_trans(3), 'wo', 'MarkerSize', 15, 'LineWidth', 2);
    text(v_trans(1)/v_trans(3), v_trans(2)/v_trans(3), ' v_{trans}', 'Color', 'w', 'FontSize', 12);
end

%% 5. Camera Calibration (K)
% Compute matrix K using orthogonality of {v_vert, v_axis, v_trans}.
% Theory: v_i' * omega * v_j = 0.
% K depends on four parameters: fx, fy, Uo, Vo (Zero skew).
% Thus omega = [w1 0 w2; 0 w3 w4; w2 w4 w5] (5 non-zero entries).
% Note: w1 != w3 because fx != fy.

fprintf('Computing Calibration Matrix K...\n');
% TODO: Solve linear system A*x = 0 for omega entries.
% K = inv(chol(omega));

%% 4. Rectification
% Find Euclidean rectification mapping H_R for a vertical plane perpendicular to axis.
% (Using K and vanishing points, or metric properties of l_inf).

fprintf('Computing Rectification...\n');
% TODO: Construct H_R
% img_rect = imwarp(img, projective2d(H_R'));

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

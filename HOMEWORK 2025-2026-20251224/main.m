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

% 2. Check for missing pieces and ask ONLY for what is missing or if re-pick is desired
choice = questdlg('How would you like to proceed?', 'Feature Extraction', 'Continue/Refine Missing', 'Start Over (Clear All)', 'Continue/Refine Missing');
if strcmp(choice, 'Start Over (Clear All)')
    lines_v = []; lines_axis = []; line_apical = []; lines_trans = [];
    arcs_A = {}; arcs_B = {};
    fprintf('Selection cleared. Starting fresh.\n');
end

% Helper function for confirmation
confirm_keep = @(name, var) ~isempty(var) && strcmp(questdlg(sprintf('%s features found. Keep them?', name), 'Review Selection', 'Keep', 'Re-pick', 'Keep'), 'Keep');

% --- Vertical Lines (Black) ---
if ~confirm_keep('Vertical', lines_v)
    fprintf('1. Select VERTICAL lines (Black in instructions).\n');
    lines_v = select_lines(img, 'Select Vertical Lines (Black). Enter to finish.', 'k');
    needs_save = true;
end

% --- Axis Parallel Lines (Green) ---
if ~confirm_keep('Axis-Parallel', lines_axis)
    fprintf('2. Select AXIS-PARALLEL lines (Green in instructions).\n');
    lines_axis = select_lines(img, 'Select Axis-Parallel Lines (Green). Enter to finish.', 'g');
    needs_save = true;
end

% --- Apical Line (Yellow) ---
if ~confirm_keep('Apical', line_apical)
    fprintf('3. Select THE APICAL line (Yellow in instructions).\n');
    temp_apical = select_lines(img, 'Select THE Apical Line (Yellow). Select ONE and Enter.', 'y');
    if ~isempty(temp_apical), line_apical = temp_apical(1); end
    needs_save = true;
end

% --- Transversal Lines (White) ---
if ~confirm_keep('Transversal', lines_trans)
    fprintf('4. Select TRANSVERSAL lines (White in instructions).\n');
    lines_trans = select_lines(img, 'Select Transversal Lines (White). Enter to finish.', 'w');
    needs_save = true;
end

% --- Arcs Family A (Cyan) ---
if ~confirm_keep('Arcs Family A', arcs_A)
    fprintf('5. Select points for DIAGONAL ARCS (Family A - Cyan).\n');
    arcs_A = select_arcs(img, 'Select Family A Arcs (Up-Right /). Enter on empty for next Step');
    needs_save = true;
end

% --- Arcs Family B (Magenta) ---
if ~confirm_keep('Arcs Family B', arcs_B)
    fprintf('6. Select points for DIAGONAL ARCS (Family B - Magenta).\n');
    arcs_B = select_arcs(img, 'Select Family B Arcs (Up-Left \). Enter twice to finish.');
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
% Incorporating the yellow apical line improves accuracy significantly.
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

% Visualize vanishing points on special figure
if ~isempty(v_vert) || ~isempty(v_axis) || ~isempty(v_trans)
    figure(2); imshow(img); title('Vanishing Points Visualization'); hold on;
    if ~isempty(v_vert), plot(v_vert(1)/v_vert(3), v_vert(2)/v_vert(3), 'ko', 'MarkerSize', 15, 'LineWidth', 2); end
    if ~isempty(v_axis), plot(v_axis(1)/v_axis(3), v_axis(2)/v_axis(3), 'go', 'MarkerSize', 15, 'LineWidth', 2); end
    if ~isempty(v_trans), plot(v_trans(1)/v_trans(3), v_trans(2)/v_trans(3), 'wo', 'MarkerSize', 15, 'LineWidth', 2); end
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

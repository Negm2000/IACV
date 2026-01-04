% IACV Homework 2025-2026
% 3D Reconstruction of a Cylindric Vault
% Updated: 2026-01-03 (Smart Loading Version)

clear; close all; clc;

%% Setup and Parameters
imageFile = 'San Maurizio.jpg';
featureFile = 'features_fixed.mat'; % Uses the repaired file

%% 1. Load Image
if exist(imageFile, 'file')
    img = imread(imageFile);
    figure(1); imshow(img); title('Input Image: San Maurizio');
    hold on;
else
    error('Image file %s not found. Please make sure it is in the current directory.', imageFile);
end

[rows, cols, ~] = size(img);

%% 2. Feature Extraction (Smart Update)
% This section checks what you HAVE vs what you NEED.

% Initialize empty variables just in case
lines_v = []; lines_axis = []; line_apical = []; lines_trans = [];
arcs_A = {}; arcs_B = {};
needs_save = false;

% 1. Try to load existing data
if exist(featureFile, 'file')
    fprintf('Loading existing features from %s...\n', featureFile);
    load(featureFile);
else
    fprintf('No feature file found. Starting fresh.\n');
end

% 2. Check for missing pieces and ask ONLY for what is missing

% --- Vertical Lines (Black) ---
if isempty(lines_v)
    fprintf('Missing VERTICAL lines.\n');
    lines_v = select_lines(img, '1. Select Vertical Lines (Black). Enter to finish.', 'k');
    needs_save = true;
end

% --- Axis Lines (Green) ---
if isempty(lines_axis)
    fprintf('Missing AXIS-PARALLEL lines.\n');
    lines_axis = select_lines(img, '2. Select Axis-Parallel Lines (Green). Enter to finish.', 'g');
    needs_save = true;
end

% --- Apical Line (Yellow) ---
if isempty(line_apical)
    fprintf('Missing APICAL line.\n');
    temp = select_lines(img, '3. Select THE Apical Line (Yellow). Enter to finish.', 'y');
    if ~isempty(temp), line_apical = temp(1); end
    needs_save = true;
end

% --- Transversal Lines (White) ---
if isempty(lines_trans)
    fprintf('Missing TRANSVERSAL lines.\n');
    lines_trans = select_lines(img, '4. Select Transversal Lines (White). Enter to finish.', 'w');
    needs_save = true;
end

% --- Arcs Family A (Up-Right /) ---
if isempty(arcs_A)
    fprintf('Missing Arcs Family A.\n');
    arcs_A = select_arcs(img, '5. Select Arcs Family A (Up-Right /). Enter to save arc, Enter again to finish.');
    needs_save = true;
end

% --- Arcs Family B (Up-Left \) ---
if isempty(arcs_B)
    fprintf('Missing Arcs Family B.\n');
    arcs_B = select_arcs(img, '6. Select Arcs Family B (Up-Left \). Enter to save arc, Enter again to finish.');
    needs_save = true;
end

% 3. Save if we added anything new
if needs_save
    save(featureFile, 'lines_v', 'lines_axis', 'line_apical', 'lines_trans', 'arcs_A', 'arcs_B');
    fprintf('All features updated and saved to %s\n', featureFile);
end

%% Visualization of All Features
figure(1); imshow(img); title('All Features Loaded'); hold on;
% Plot lines
if ~isempty(lines_v), for i=1:length(lines_v), plot([lines_v(i).p1(1) lines_v(i).p2(1)], [lines_v(i).p1(2) lines_v(i).p2(2)], 'k-','LineWidth',2); end; end
if ~isempty(lines_axis), for i=1:length(lines_axis), plot([lines_axis(i).p1(1) lines_axis(i).p2(1)], [lines_axis(i).p1(2) lines_axis(i).p2(2)], 'g-','LineWidth',2); end; end
if ~isempty(line_apical), plot([line_apical.p1(1) line_apical.p2(1)], [line_apical.p1(2) line_apical.p2(2)], 'y-','LineWidth',3); end
if ~isempty(lines_trans), for i=1:length(lines_trans), plot([lines_trans(i).p1(1) lines_trans(i).p2(1)], [lines_trans(i).p1(2) lines_trans(i).p2(2)], 'w-','LineWidth',2); end; end
% Plot Arcs
for i=1:length(arcs_A)
    pts=arcs_A{i};
    plot(pts(:,1), pts(:,2), 'c.-','LineWidth',1.5, 'MarkerSize', 10);
end % Cyan for A
for i=1:length(arcs_B)
    pts=arcs_B{i};
    plot(pts(:,1), pts(:,2), 'm.-','LineWidth',1.5, 'MarkerSize', 10);
end % Magenta for B

%% 3. Vanishing Points
fprintf('Computing Vanishing Points...\n');
get_line = @(p1, p2) cross([p1, 1], [p2, 1])';

% v_vert (Vertical)
lines_hom = {};
if ~isempty(lines_v)
    for i=1:length(lines_v), lines_hom{end+1} = get_line(lines_v(i).p1, lines_v(i).p2); end
    v_vert = compute_vanishing_point(lines_hom);
else
    v_vert = [];
end

% v_trans (Transversal)
lines_hom = {};
if ~isempty(lines_trans)
    for i=1:length(lines_trans), lines_hom{end+1} = get_line(lines_trans(i).p1, lines_trans(i).p2); end
    v_trans = compute_vanishing_point(lines_hom);
else
    v_trans = [];
end

% v_axis (IMPROVED: Green + Yellow Combined)
lines_hom = {};
has_axis_info = false;
if ~isempty(lines_axis)
    for i=1:length(lines_axis), lines_hom{end+1} = get_line(lines_axis(i).p1, lines_axis(i).p2); end
    has_axis_info = true;
end
if ~isempty(line_apical)
    lines_hom{end+1} = get_line(line_apical.p1, line_apical.p2);
    has_axis_info = true;
end

if has_axis_info
    v_axis = compute_vanishing_point(lines_hom);
else
    v_axis = [];
end

if ~isempty(v_vert),  fprintf('v_vert  = [%.4f, %.4f, %.4f]\n', v_vert); end
if ~isempty(v_axis),  fprintf('v_axis  = [%.4f, %.4f, %.4f]\n', v_axis); end
if ~isempty(v_trans), fprintf('v_trans = [%.4f, %.4f, %.4f]\n', v_trans); end

% Compute Vanishing Line of planes perpendicular to the cylinder axis
if ~isempty(v_vert) && ~isempty(v_trans)
    l_inf_perp = cross(v_vert, v_trans);
    l_inf_perp = l_inf_perp / norm(l_inf_perp(1:2));
    fprintf('Vanishing line l_inf_perp = [%.4f, %.4f, %.4f]\n', l_inf_perp);
else
    l_inf_perp = [];
end

%% 5. Camera Calibration (K)
% Compute matrix K using orthogonality of {v_vert, v_axis, v_trans}.
% Theory: v_i' * omega * v_j = 0.
% K depends on four parameters: fx, fy, Uo, Vo (Zero skew).
% Thus omega = [w1 0 w2; 0 w3 w4; w2 w4 w5] (5 non-zero entries).
% Note: w1 != w3 because fx != fy.

fprintf('Computing Calibration Matrix K...\n');
% TODO: Implement solver for omega and extract K

%% 4. Rectification
% Find Euclidean rectification mapping H_R for a vertical plane perpendicular to axis.
fprintf('Computing Rectification...\n');
% TODO: Construct H_R

%% 6. 3D Reconstruction
% BACK-PROJECTION AND INTERSECTION WITH CYLINDER
fprintf('Performing 3D Reconstruction...\n');
% TODO: Implementation of ray-cylinder intersection logic.

%% Part 2 - Visualization

%% 7. Plotting
figure('Name', '3D Reconstruction');
% TODO: plot3(...)
grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');

disp('Processing Complete.');

%% Helper Functions

function vp = compute_vanishing_point(lines_hom)
if length(lines_hom) < 2
    warning('Not enough lines for VP. Need at least 2.');
    vp = [0; 0; 0];
    return;
end
A = zeros(length(lines_hom), 3);
for k = 1:length(lines_hom), A(k, :) = lines_hom{k}'; end
[~, ~, V] = svd(A);
vp = V(:, end);
if abs(vp(3)) > 1e-9, vp = vp / vp(3); end
vp = vp';
end

function lines = select_lines(img, promptStr, colorCode)
% Robust line selection with Zoom support and single-struct storage.
hFig = figure('Name', promptStr); imshow(img); hold on;
title({promptStr, 'MOUSE: Zoom/Pan with Toolbar tools', 'PICK: Toggle tool OFF to click', '[u]: Undo, [Enter]: Finish'});

lines = struct('p1', {}, 'p2', {});
plotHandles = {};
tempHandle = [];

while true
    k = waitforbuttonpress;
    if k == 0 % Mouse Button
        % Check if zoom/pan active
        z = zoom(hFig); p = pan(hFig);
        if strcmp(z.Enable, 'on') || strcmp(p.Enable, 'on'), continue; end

        currPt = get(gca, 'CurrentPoint');
        pt = currPt(1,1:2);

        if isempty(tempHandle)
            % High precision crosshair for point 1
            tempHandle = plot(pt(1), pt(2), [colorCode '+'], 'MarkerSize', 12, 'LineWidth', 1.5);
            p1_val = pt;
        else
            % Point 2
            p2_val = pt;
            % Store in single struct entry
            idx = length(lines) + 1;
            lines(idx).p1 = p1_val;
            lines(idx).p2 = p2_val;

            delete(tempHandle); tempHandle = [];
            % visual feedback
            hL = plot([p1_val(1) p2_val(1)], [p1_val(2) p2_val(2)], [colorCode '-'], 'LineWidth', 2);
            hP1 = plot(p1_val(1), p1_val(2), [colorCode '+'], 'MarkerSize', 8);
            hP2 = plot(p2_val(1), p2_val(2), [colorCode '+'], 'MarkerSize', 8);
            plotHandles{end+1} = [hL, hP1, hP2];
            fprintf('Line added.\n');
        end

    elseif k == 1 % Keyboard
        key = get(hFig, 'CurrentCharacter');
        if isempty(key), continue; end

        if key == 13 % Enter
            if isempty(tempHandle), break; else fprintf('Finish current line first.\n'); end
        elseif key == 'u' % Undo
            if ~isempty(tempHandle)
                delete(tempHandle); tempHandle = [];
                fprintf('Point 1 cancelled.\n');
            elseif ~isempty(lines)
                delete(plotHandles{end});
                plotHandles(end) = [];
                lines(end) = [];
                fprintf('Last line undone.\n');
            end
        end
    end
end
close(hFig);
end

function arcs = select_arcs(img, promptStr)
% Arc selection logic (Red points for Family A/B)
hFig = figure('Name', promptStr); imshow(img); hold on;
title({promptStr, 'MOUSE: Click to pick (Ensure Tool is OFF)', '[u/Backspace]: Undo, [Enter]: Save Arc/Finish'});

arcs = {};
currentArcPts = [];
currentArcPlot = [];

while true
    k = waitforbuttonpress;
    if k == 0 % Click
        z = zoom(hFig); p = pan(hFig);
        if strcmp(z.Enable, 'on') || strcmp(p.Enable, 'on'), continue; end

        currPt = get(gca, 'CurrentPoint');
        pt = currPt(1,1:2);
        currentArcPts = [currentArcPts; pt];

        if ~isempty(currentArcPlot), delete(currentArcPlot); end
        currentArcPlot = plot(currentArcPts(:,1), currentArcPts(:,2), 'r.-', 'MarkerSize', 12, 'LineWidth', 1.5);

    elseif k == 1 % Key
        key = get(hFig, 'CurrentCharacter');
        if isempty(key), continue; end

        if key == 13 % Enter
            if isempty(currentArcPts), break; end
            arcs{end+1} = currentArcPts;
            currentArcPts = []; currentArcPlot = [];
            fprintf('Arc saved.\n');
        elseif key == 'u' || key == 8 % 'u' or Backspace
            if ~isempty(currentArcPts)
                currentArcPts(end,:) = [];
                if ~isempty(currentArcPlot), delete(currentArcPlot); end
                if ~isempty(currentArcPts)
                    currentArcPlot = plot(currentArcPts(:,1), currentArcPts(:,2), 'r.-', 'MarkerSize', 12);
                else
                    currentArcPlot = [];
                end
            elseif ~isempty(arcs)
                % Undo whole arc not implemented in user snippet but good for consistency
                arcs(end) = [];
                fprintf('Last arc removed.\n');
            end
        end
    end
end
close(hFig);
end

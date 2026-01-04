% IACV Homework 2025-2026
% 3D Reconstruction of a Cylindric Vault
% Updated: 2026-01-03 (Smart Loading Version)

clear; close all; clc;

%% Setup and Parameters
imageFile = 'San Maurizio.jpg';
featureFile = 'features_fixed.mat'; % Uses your repaired file

%% 1. Load Image
if exist(imageFile, 'file')
    img = imread(imageFile);
    figure(1); imshow(img); title('Input Image: San Maurizio');
    hold on;
else
    error('Image file %s not found.', imageFile);
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
% [THIS WAS MISSING BEFORE]
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
for i=1:length(arcs_A), pts=arcs_A{i}; plot(pts(:,1), pts(:,2), 'c.-','LineWidth',1.5); end % Cyan for A
for i=1:length(arcs_B), pts=arcs_B{i}; plot(pts(:,1), pts(:,2), 'm.-','LineWidth',1.5); end % Magenta for B

%% 3. Vanishing Points
fprintf('Computing Vanishing Points...\n');
get_line = @(p1, p2) cross([p1, 1], [p2, 1])';

% v_vert (Vertical)
lines_hom = {};
for i=1:length(lines_v), lines_hom{end+1} = get_line(lines_v(i).p1, lines_v(i).p2); end
v_vert = compute_vanishing_point(lines_hom);

% v_trans (Transversal)
lines_hom = {};
for i=1:length(lines_trans), lines_hom{end+1} = get_line(lines_trans(i).p1, lines_trans(i).p2); end
v_trans = compute_vanishing_point(lines_hom);

% v_axis (IMPROVED: Green + Yellow Combined)
lines_hom = {};
for i=1:length(lines_axis), lines_hom{end+1} = get_line(lines_axis(i).p1, lines_axis(i).p2); end
if ~isempty(line_apical), lines_hom{end+1} = get_line(line_apical.p1, line_apical.p2); end
v_axis = compute_vanishing_point(lines_hom);

fprintf('v_vert  = [%.4f, %.4f, %.4f]\n', v_vert);
fprintf('v_axis  = [%.4f, %.4f, %.4f]\n', v_axis);
fprintf('v_trans = [%.4f, %.4f, %.4f]\n', v_trans);

%% 4 & 5 & 6 (Placeholders)
% [Keep your TODOs here]

%% Helper Functions

function vp = compute_vanishing_point(lines_hom)
    if length(lines_hom) < 2, warning('Not enough lines for VP'); vp=[0;0;0]; return; end
    A = zeros(length(lines_hom), 3);
    for k = 1:length(lines_hom), A(k, :) = lines_hom{k}'; end
    [~, ~, V] = svd(A);
    vp = V(:, end);
    if abs(vp(3)) > 1e-9, vp = vp / vp(3); end
    vp = vp';
end

function lines = select_lines(img, promptStr, colorCode)
% [Use your FIXED version of select_lines here]
% Make sure to keep the "lines(idx).p1 = p1..." fix you added!
% I am including the structure wrapper just to be safe:
hFig = figure('Name', promptStr); imshow(img); hold on;
title(promptStr);
lines = struct('p1', {}, 'p2', {});
tempHandle = [];
while true
    k = waitforbuttonpress;
    if k == 1 % Key
        key = get(hFig, 'CurrentCharacter');
        if key==13, break; end % Enter
        if key=='u' && ~isempty(lines) % Undo
           delete(findobj(gca, 'Color', colorCode)); % Clear visual
           lines(end) = []; % Remove data
           % Redraw remaining lines
           for i=1:length(lines)
               plot([lines(i).p1(1) lines(i).p2(1)], [lines(i).p1(2) lines(i).p2(2)], [colorCode '-'], 'LineWidth', 2);
           end
        end
    elseif k == 0 % Click
        currPt = get(gca, 'CurrentPoint'); p = currPt(1,1:2);
        if isempty(tempHandle)
            tempHandle = plot(p(1), p(2), [colorCode '+']);
            p1 = p;
        else
            delete(tempHandle); tempHandle = [];
            p2 = p;
            lines(end+1).p1 = p1; lines(end).p2 = p2; % The Safe Way
            plot([p1(1) p2(1)], [p1(2) p2(2)], [colorCode '-'], 'LineWidth', 2);
        end
    end
end
close(hFig);
end

function arcs = select_arcs(img, promptStr)
% [Use your existing select_arcs, it was fine]
hFig = figure('Name', promptStr); imshow(img); hold on;
title(promptStr);
arcs = {}; currentArcPts = []; currentArcPlot = [];
while true
    k = waitforbuttonpress;
    if k == 0 % Click
         currPt = get(gca, 'CurrentPoint'); p = currPt(1,1:2);
         currentArcPts = [currentArcPts; p];
         if ~isempty(currentArcPlot), delete(currentArcPlot); end
         currentArcPlot = plot(currentArcPts(:,1), currentArcPts(:,2), 'r.-');
    elseif k == 1
        key = get(hFig, 'CurrentCharacter');
        if key == 13
            if isempty(currentArcPts), break; end
            arcs{end+1} = currentArcPts;
            currentArcPts = []; currentArcPlot = [];
        end
    end
end
close(hFig);
end
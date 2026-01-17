% feature_picker.m - Interactive Feature Selection for IACV Homework
% This script handles all feature picking and saves to a .mat file.
% Run this separately to pick/re-pick features, then run main.m for analysis.

clear; close all; clc;

%% Configuration
imageFile = 'San Maurizio.jpg';
featureFile = 'features_fixed.mat';

%% Load Image
if exist(imageFile, 'file')
    img = imread(imageFile);
else
    error('Image file %s not found.', imageFile);
end

%% Initialize feature variables
lines_v = [];
lines_axis = [];
line_apical = [];
lines_trans = [];
arcs_A = {};
arcs_B = {};

%% Load existing features if available
if exist(featureFile, 'file')
    fprintf('Loading existing features from %s...\n', featureFile);
    load(featureFile);
end

%% Ask which features to pick/re-pick
repickList = [];
prompt = {['Select features to PICK/RE-PICK (type numbers e.g. "1 3 5" or leave empty to KEEP ALL):', newline, ...
    '1. Vertical (Black)', newline, ...
    '2. Axis-Parallel (Green)', newline, ...
    '3. Apical (Yellow)', newline, ...
    '4. Transversal (White)', newline, ...
    '5. Arcs Family A (Cyan)', newline, ...
    '6. Arcs Family B (Magenta)']};
answer = inputdlg(prompt, 'Feature Picker', [1 100]);

if ~isempty(answer) && ~isempty(answer{1})
    repickList = str2num(answer{1});
end

%% Pick features as needed
if isempty(lines_v) || ismember(1, repickList)
    fprintf('1. Select VERTICAL lines (Black).\n');
    lines_v = select_lines(img, 'Select Vertical Lines (Black). Enter to finish.', 'k');
end

if isempty(lines_axis) || ismember(2, repickList)
    fprintf('2. Select AXIS-PARALLEL lines (Green).\n');
    lines_axis = select_lines(img, 'Select Axis-Parallel Lines (Green). Enter to finish.', 'g');
end

if isempty(line_apical) || ismember(3, repickList)
    fprintf('3. Select THE APICAL line (Yellow).\n');
    temp_apical = select_lines(img, 'Select THE Apical Line (Yellow). Select ONE and Enter.', 'y');
    if ~isempty(temp_apical), line_apical = temp_apical(1); end
end

if isempty(lines_trans) || ismember(4, repickList)
    fprintf('4. Select TRANSVERSAL lines (White).\n');
    lines_trans = select_lines(img, 'Select Transversal Lines (White). Enter to finish.', 'w');
end

if isempty(arcs_A) || ismember(5, repickList)
    fprintf('5. Select points for DIAGONAL ARCS (Family A - Cyan).\n');
    arcs_A = select_arcs(img, 'Select Family A Arcs (Up-Right /). Enter on empty to Finish');
end

if isempty(arcs_B) || ismember(6, repickList)
    fprintf('6. Select points for DIAGONAL ARCS (Family B - Magenta).\n');
    arcs_B = select_arcs(img, 'Select Family B Arcs (Up-Left \). Enter twice to Finish.');
end

%% Save features
save(featureFile, 'lines_v', 'lines_axis', 'line_apical', 'lines_trans', 'arcs_A', 'arcs_B');
fprintf('Features saved to %s\n', featureFile);

%% Display all features
figure(1); imshow(img); title('All Features'); hold on;
if ~isempty(lines_v), for i=1:length(lines_v), plot([lines_v(i).p1(1) lines_v(i).p2(1)], [lines_v(i).p1(2) lines_v(i).p2(2)], 'k-', 'LineWidth', 2); end; end
if ~isempty(lines_axis), for i=1:length(lines_axis), plot([lines_axis(i).p1(1) lines_axis(i).p2(1)], [lines_axis(i).p1(2) lines_axis(i).p2(2)], 'g-', 'LineWidth', 2); end; end
if ~isempty(line_apical), plot([line_apical.p1(1) line_apical.p2(1)], [line_apical.p1(2) line_apical.p2(2)], 'y-', 'LineWidth', 3); end
if ~isempty(lines_trans), for i=1:length(lines_trans), plot([lines_trans(i).p1(1) lines_trans(i).p2(1)], [lines_trans(i).p1(2) lines_trans(i).p2(2)], 'w-', 'LineWidth', 2); end; end
for i=1:length(arcs_A), pts=arcs_A{i}; plot(pts(:,1), pts(:,2), 'c.-', 'MarkerSize', 10); end
for i=1:length(arcs_B), pts=arcs_B{i}; plot(pts(:,1), pts(:,2), 'm.-', 'MarkerSize', 10); end

fprintf('\n=== Feature Summary ===\n');
fprintf('Vertical lines: %d\n', length(lines_v));
fprintf('Axis-parallel lines: %d\n', length(lines_axis));
fprintf('Apical line: %s\n', ternary(~isempty(line_apical), 'Yes', 'No'));
fprintf('Transversal lines: %d\n', length(lines_trans));
fprintf('Arc family A: %d arcs\n', length(arcs_A));
fprintf('Arc family B: %d arcs\n', length(arcs_B));
fprintf('\nRun main.m to perform the analysis.\n');

%% Helper function
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

%% Line selection function
function lines = select_lines(img, promptStr, colorCode)
if nargin < 3, colorCode = 'y'; end
hFig = figure('Name', promptStr); imshow(img); hold on;
title({promptStr, 'MOUSE: Zoom/Pan with Toolbar tools', ...
    'ACTION: Toggle tools OFF to click/pick points', ...
    '[u]: Undo, [Enter]: Finish category'});

lines = struct('p1', {}, 'p2', {});
plotHandles = {};
tempHandle = [];

while true
    k = waitforbuttonpress;

    if k == 0
        zoomObj = zoom(hFig);
        panObj = pan(hFig);
        if strcmp(zoomObj.Enable, 'on') || strcmp(panObj.Enable, 'on')
            continue;
        end

        currPt = get(gca, 'CurrentPoint');
        x = currPt(1,1);
        y = currPt(1,2);

        if isempty(tempHandle)
            tempHandle = plot(x, y, [colorCode '+'], 'MarkerSize', 12, 'LineWidth', 1.5);
            p1 = [x, y];
            fprintf('Point 1 set at (%.2f, %.2f)\n', x, y);
        else
            p2 = [x, y];
            idx = length(lines) + 1;
            lines(idx).p1 = p1;
            lines(idx).p2 = p2;
            delete(tempHandle);
            tempHandle = [];
            hLine = plot([p1(1), p2(1)], [p1(2), p2(2)], [colorCode '-'], 'LineWidth', 2);
            hP1 = plot(p1(1), p1(2), [colorCode '+'], 'MarkerSize', 8);
            hP2 = plot(p2(1), p2(2), [colorCode '+'], 'MarkerSize', 8);
            plotHandles{end+1} = [hLine, hP1, hP2];
            fprintf('Line added.\n');
        end

    elseif k == 1
        key = get(hFig, 'CurrentCharacter');
        if isempty(key), continue; end

        if key == 13
            if isempty(tempHandle)
                break;
            else
                fprintf('Current line incomplete. Pick second point or [u] to cancel.\n');
            end
        elseif key == 'u'
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

%% Arc selection function
function arcs = select_arcs(img, promptStr)
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
    if k == 0
        zoomObj = zoom(hFig); panObj = pan(hFig);
        if strcmp(zoomObj.Enable, 'on') || strcmp(panObj.Enable, 'on'), continue; end

        currPt = get(gca, 'CurrentPoint');
        x = currPt(1,1); y = currPt(1,2);
        currentArcPts = [currentArcPts; x, y];
        if ~isempty(currentArcPlot), delete(currentArcPlot); end
        currentArcPlot = plot(currentArcPts(:,1), currentArcPts(:,2), 'r.-', 'MarkerSize', 12, 'LineWidth', 1.5);
        fprintf('Point picked at (%.1f, %.1f)\n', x, y);

    elseif k == 1
        key = get(hFig, 'CurrentCharacter');
        if isempty(key), continue; end

        if key == 13
            if isempty(currentArcPts)
                break;
            else
                arcs{end+1} = currentArcPts;
                arcPlots{end+1} = currentArcPlot;
                currentArcPts = [];
                currentArcPlot = [];
                fprintf('Arc saved.\n');
            end
        elseif key == 'u' || key == 8
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

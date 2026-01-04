% NUCLEAR_REPAIR.m
% Fixes "Index exceeds array bounds" by removing orphaned points caused by Undo.

clear; clc;
filename = 'features.mat'; % The broken file

if ~exist(filename, 'file')
    error('File features.mat not found.');
end
load(filename);

% --- The Smart Stitching Function ---
function fixed = smart_stitch(broken)
    fixed = struct('p1', {}, 'p2', {});
    if isempty(broken), return; end
    
    i = 1;
    while i < length(broken)
        curr = broken(i);
        next = broken(i+1);
        
        % Check if Current is a Start (p1) and Next is an End (p2)
        has_p1 = ~isempty(curr.p1);
        next_has_p2 = ~isempty(next.p2);
        
        if has_p1 && next_has_p2
            % MATCH FOUND: It's a valid pair. Save it.
            fixed(end+1).p1 = curr.p1;
            fixed(end).p2 = next.p2;
            i = i + 2; % Skip both
        else
            % MISMATCH: This means 'curr' was an orphan (Undo button residue).
            % Skip it and try to match the next one.
            i = i + 1;
        end
    end
end
% ------------------------------------

fprintf('Repairing data with smart logic...\n');

% Repair all categories
if exist('lines_v', 'var'), lines_v = smart_stitch(lines_v); end
if exist('lines_axis', 'var'), lines_axis = smart_stitch(lines_axis); end
if exist('line_apical', 'var'), line_apical = smart_stitch(line_apical); end
if exist('lines_trans', 'var'), lines_trans = smart_stitch(lines_trans); end

% Save to a CLEAN file
save('features_fixed.mat', 'lines_v', 'lines_axis', 'line_apical', 'lines_trans', 'arcs_A');
fprintf('SUCCESS. Fixed data saved to "features_fixed.mat".\n');
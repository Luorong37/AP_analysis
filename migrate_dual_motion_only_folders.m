% migrate_dual_motion_only_folders
%
% Recursively migrate legacy Dual_analysis3_motion_only result folders into
% the normal Dual_analysis3 result tree. This treats motion_only as a stop
% point inside Dual_analysis3, not as a parallel analysis module.
%
% Usage:
%   root_path = 'E:\path\to\Methods_or_Rec_folder';
%   run('migrate_dual_motion_only_folders.m');
%
% Optional:
%   dry_run = true;   % preview actions without moving folders/files

if ~exist('root_path', 'var') || isempty(root_path)
    root_path = pwd;
end
if ~exist('dry_run', 'var') || isempty(dry_run)
    dry_run = false;
end

legacy_folder_name = 'Dual_analysis3_motion_only';
target_folder_name = 'Dual_analysis3';

if ~isfolder(root_path)
    error('root_path does not exist or is not a folder: %s', root_path);
end

fprintf('\n============================================================\n');
fprintf('Migrate Dual_analysis3 Motion-Only Folders\n');
fprintf('============================================================\n');
fprintf('Root: %s\n', root_path);
fprintf('Legacy folder: %s\n', legacy_folder_name);
fprintf('Target folder: %s\n', target_folder_name);
fprintf('Dry run: %d\n', logical(dry_run));

legacy_dirs = find_named_folders_recursive(root_path, legacy_folder_name);
fprintf('Found legacy folders: %d\n', numel(legacy_dirs));

migration_log = repmat(struct( ...
    'legacy_path', "", ...
    'target_path', "", ...
    'status', "", ...
    'moved_items', 0, ...
    'renamed_conflicts', 0, ...
    'message', ""), 0, 1);

for idx = 1:numel(legacy_dirs)
    legacy_path = char(legacy_dirs(idx));
    parent_path = fileparts(legacy_path);
    target_path = fullfile(parent_path, target_folder_name);

    fprintf('\n[%d/%d]\n', idx, numel(legacy_dirs));
    fprintf('  From: %s\n', legacy_path);
    fprintf('  To:   %s\n', target_path);

    try
        if ~isfolder(target_path)
            if dry_run
                fprintf('  Action: would rename folder.\n');
            else
                movefile(legacy_path, target_path);
                fprintf('  Action: renamed folder.\n');
            end
            migration_log(end+1, 1) = struct( ...
                'legacy_path', string(legacy_path), ...
                'target_path', string(target_path), ...
                'status', "renamed_folder", ...
                'moved_items', 1, ...
                'renamed_conflicts', 0, ...
                'message', "");
        else
            [moved_items, renamed_conflicts] = merge_folder_contents(legacy_path, target_path, dry_run);
            if dry_run
                fprintf('  Action: would merge %d item(s), conflicts renamed: %d.\n', moved_items, renamed_conflicts);
            else
                remove_folder_if_empty(legacy_path);
                fprintf('  Action: merged %d item(s), conflicts renamed: %d.\n', moved_items, renamed_conflicts);
            end
            migration_log(end+1, 1) = struct( ...
                'legacy_path', string(legacy_path), ...
                'target_path', string(target_path), ...
                'status', "merged_into_existing_folder", ...
                'moved_items', moved_items, ...
                'renamed_conflicts', renamed_conflicts, ...
                'message', "");
        end
    catch ME
        migration_log(end+1, 1) = struct( ...
            'legacy_path', string(legacy_path), ...
            'target_path', string(target_path), ...
            'status', "failed", ...
            'moved_items', 0, ...
            'renamed_conflicts', 0, ...
            'message', string(getReport(ME, 'extended', 'hyperlinks', 'off')));
        fprintf(2, '  Failed: %s\n', ME.message);
    end
end

fprintf('\nMigration complete.\n');
fprintf('Completed entries: %d\n', sum(string({migration_log.status}) ~= "failed"));
fprintf('Failed entries: %d\n', sum(string({migration_log.status}) == "failed"));
fprintf('The variable migration_log is available in the workspace.\n');

function folder_paths = find_named_folders_recursive(root_path, folder_name)
folder_paths = strings(0, 1);
items = dir(root_path);
items = items([items.isdir]);
items = items(~ismember({items.name}, {'.', '..'}));
for item_idx = 1:numel(items)
    current_path = fullfile(items(item_idx).folder, items(item_idx).name);
    if strcmp(items(item_idx).name, folder_name)
        folder_paths(end+1, 1) = string(current_path);
    else
        child_paths = find_named_folders_recursive(current_path, folder_name);
        folder_paths = [folder_paths; child_paths]; %#ok<AGROW>
    end
end
end

function [moved_items, renamed_conflicts] = merge_folder_contents(source_folder, target_folder, dry_run)
moved_items = 0;
renamed_conflicts = 0;

if ~isfolder(target_folder) && ~dry_run
    mkdir(target_folder);
end

items = dir(source_folder);
items = items(~ismember({items.name}, {'.', '..'}));
for item_idx = 1:numel(items)
    source_path = fullfile(items(item_idx).folder, items(item_idx).name);
    destination_path = fullfile(target_folder, items(item_idx).name);
    if exist(destination_path, 'file') || isfolder(destination_path)
        destination_path = make_unique_destination(target_folder, items(item_idx).name, items(item_idx).isdir);
        renamed_conflicts = renamed_conflicts + 1;
    end

    moved_items = moved_items + 1;
    fprintf('    %s -> %s\n', source_path, destination_path);
    if ~dry_run
        movefile(source_path, destination_path);
    end
end
end

function destination_path = make_unique_destination(target_folder, original_name, is_directory)
if is_directory
    base_name = original_name;
    ext_part = '';
else
    [~, name_part, ext_part] = fileparts(original_name);
    base_name = name_part;
end

time_tag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
candidate_name = sprintf('%s_migrated_%s%s', base_name, time_tag, ext_part);
destination_path = fullfile(target_folder, candidate_name);
counter = 1;
while exist(destination_path, 'file') || isfolder(destination_path)
    candidate_name = sprintf('%s_migrated_%s_%03d%s', base_name, time_tag, counter, ext_part);
    destination_path = fullfile(target_folder, candidate_name);
    counter = counter + 1;
end
end

function remove_folder_if_empty(folder_path)
items = dir(folder_path);
items = items(~ismember({items.name}, {'.', '..'}));
if isempty(items)
    rmdir(folder_path);
end
end

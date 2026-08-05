function entries = discover_record_cycles(scope_path, cycle_filter, allow_nested_records)
%DISCOVER_RECORD_CYCLES Discover deterministic Record/Cycle jobs.
%   ENTRY fields: label, cycle_name, record_name, record_path, cycle_path.

scope_path = char(string(scope_path));
if ~isfolder(scope_path)
    error('AAA:Record:MissingScope', 'Record scope does not exist: %s', scope_path);
end
if nargin < 2 || isempty(cycle_filter)
    cycle_filter = strings(0, 1);
end
if nargin < 3 || isempty(allow_nested_records)
    allow_nested_records = false;
end
cycle_filter = string(cycle_filter);
cycle_filter = cycle_filter(:);

record_paths = string(scope_path);
direct_cycles = cycle_directories(scope_path);
if isempty(direct_cycles) && logical(allow_nested_records)
    rec_dirs = dir(fullfile(scope_path, 'Rec*'));
    rec_dirs = rec_dirs([rec_dirs.isdir]);
    record_paths = strings(0, 1);
    for idx = 1:numel(rec_dirs)
        candidate = string(fullfile(rec_dirs(idx).folder, rec_dirs(idx).name));
        if ~isempty(cycle_directories(candidate))
            record_paths(end+1, 1) = candidate; %#ok<AGROW>
        end
    end
end

template = struct( ...
    'label', "", ...
    'cycle_name', "", ...
    'record_name', "", ...
    'record_path', "", ...
    'cycle_path', "");
entries = repmat(template, 0, 1);
for record_idx = 1:numel(record_paths)
    record_path = char(record_paths(record_idx));
    [~, record_name] = fileparts(record_path);
    cycle_dirs = cycle_directories(record_path);
    for cycle_idx = 1:numel(cycle_dirs)
        cycle_name = string(cycle_dirs(cycle_idx).name);
        if numel(record_paths) == 1
            label = cycle_name;
        else
            label = string(record_name) + "/" + cycle_name;
        end
        if ~isempty(cycle_filter) ...
                && ~ismember(cycle_name, cycle_filter) ...
                && ~ismember(label, cycle_filter)
            continue;
        end
        entry = template;
        entry.label = label;
        entry.cycle_name = cycle_name;
        entry.record_name = string(record_name);
        entry.record_path = string(record_path);
        entry.cycle_path = string(fullfile(cycle_dirs(cycle_idx).folder, cycle_dirs(cycle_idx).name));
        entries(end+1, 1) = entry; %#ok<AGROW>
    end
end

if isempty(entries)
    error('AAA:Record:NoCycles', ...
        'No Cycle* folders remain under the requested scope and filter: %s', scope_path);
end
end

function cycle_dirs = cycle_directories(record_path)
cycle_dirs = dir(fullfile(char(record_path), 'Cycle*'));
cycle_dirs = cycle_dirs([cycle_dirs.isdir]);
if isempty(cycle_dirs)
    return;
end
names = string({cycle_dirs.name});
numbers = nan(size(names));
for idx = 1:numel(names)
    token = regexp(names(idx), '^Cycle(\d+)', 'tokens', 'once');
    if ~isempty(token)
        numbers(idx) = str2double(token{1});
    end
end
numbers(isnan(numbers)) = Inf;
sort_table = table(numbers(:), lower(names(:)), (1:numel(names))', ...
    'VariableNames', {'number','name','original_index'});
sort_table = sortrows(sort_table, {'number','name','original_index'});
cycle_dirs = cycle_dirs(sort_table.original_index);
end

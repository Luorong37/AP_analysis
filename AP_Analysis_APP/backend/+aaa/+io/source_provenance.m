function source = source_provenance()
%SOURCE_PROVENANCE Describe the exact AAA source used for a result manifest.
% This function intentionally does not launch Git. The repository commit is
% read directly from .git/HEAD, while a deterministic SHA-256 source-bundle
% digest identifies the production files even when the worktree is modified.

app_root = locate_app_root();
warnings = strings(0,1);

source = struct( ...
    'version', read_version_file(app_root), ...
    'app_version', read_version_file(app_root), ...
    'git_commit', "", ...
    'git_branch', "", ...
    'git_worktree_status', "not_checked", ...
    'repository_root', "", ...
    'app_root', string(app_root), ...
    'matlab_version', string(version), ...
    'matlab_release', string(version('-release')), ...
    'dependencies', struct('NoRMCorre', "0.1.1"), ...
    'hashes', empty_hash_record(), ...
    'captured_at', datetime('now'), ...
    'warnings', strings(0,1));

try
    [source.git_commit, source.git_branch, source.repository_root, git_warning] = ...
        read_git_metadata(app_root);
    warnings = [warnings; git_warning(:)];
catch ME
    warnings(end+1,1) = "Git metadata unavailable: " + string(ME.message);
end

try
    [bundle_hash, included_files] = hash_production_source(app_root);
    source.hashes = struct( ...
        'algorithm', "SHA-256", ...
        'source_bundle', bundle_hash, ...
        'file_count', numel(included_files), ...
        'included_files', included_files);
catch ME
    warnings(end+1,1) = "Source hash unavailable: " + string(ME.message);
end

source.warnings = warnings(strlength(warnings) > 0);
end

function app_root = locate_app_root()
io_dir = fileparts(mfilename('fullpath'));
aaa_dir = fileparts(io_dir);
backend_dir = fileparts(aaa_dir);
app_root = fileparts(backend_dir);
end

function app_version = read_version_file(app_root)
version_file = fullfile(app_root, 'VERSION');
if ~isfile(version_file)
    app_version = "unknown";
    return;
end
app_version = strip(string(fileread(version_file)));
if strlength(app_version) == 0
    app_version = "unknown";
end
end

function value = empty_hash_record()
value = struct( ...
    'algorithm', "SHA-256", ...
    'source_bundle', "", ...
    'file_count', 0, ...
    'included_files', strings(0,1));
end

function [commit, branch, repository_root, warnings] = read_git_metadata(start_dir)
commit = "";
branch = "";
repository_root = "";
warnings = strings(0,1);

[git_dir, repo_root] = find_git_directory(start_dir);
if strlength(git_dir) == 0
    warnings(end+1,1) = "No .git metadata was found above the App folder.";
    return;
end
repository_root = repo_root;

head_file = fullfile(git_dir, 'HEAD');
if ~isfile(head_file)
    warnings(end+1,1) = "The repository HEAD file is missing.";
    return;
end
head = strip(string(fileread(head_file)));
if startsWith(head, "ref:")
    ref_name = strip(extractAfter(head, "ref:"));
    branch = erase(ref_name, "refs/heads/");
    ref_file = fullfile(git_dir, char(replace(ref_name, '/', filesep)));
    if isfile(ref_file)
        commit = strip(string(fileread(ref_file)));
    else
        commit = read_packed_ref(git_dir, ref_name);
    end
else
    branch = "detached";
    commit = head;
end

if strlength(commit) == 0
    warnings(end+1,1) = "The current Git commit could not be resolved.";
end
end

function [git_dir, repository_root] = find_git_directory(start_dir)
git_dir = "";
repository_root = "";
cursor = char(start_dir);
while true
    marker = fullfile(cursor, '.git');
    if isfolder(marker)
        git_dir = string(marker);
        repository_root = string(cursor);
        return;
    end
    if isfile(marker)
        marker_text = strip(string(fileread(marker)));
        if startsWith(marker_text, "gitdir:")
            candidate = strip(extractAfter(marker_text, "gitdir:"));
            if ~is_absolute_path(candidate)
                candidate = string(fullfile(cursor, char(candidate)));
            end
            git_dir = string(candidate);
            repository_root = string(cursor);
            return;
        end
    end
    parent = fileparts(cursor);
    if strcmp(parent, cursor)
        return;
    end
    cursor = parent;
end
end

function tf = is_absolute_path(value)
value = char(value);
tf = ~isempty(regexp(value, '^[A-Za-z]:[\\/]', 'once')) || ...
    startsWith(value, '\\') || startsWith(value, '/');
end

function commit = read_packed_ref(git_dir, ref_name)
commit = "";
packed_file = fullfile(git_dir, 'packed-refs');
if ~isfile(packed_file)
    return;
end
lines = splitlines(string(fileread(packed_file)));
for idx = 1:numel(lines)
    line = strip(lines(idx));
    if strlength(line) == 0 || startsWith(line, "#") || startsWith(line, "^")
        continue;
    end
    parts = split(line);
    if numel(parts) >= 2 && parts(2) == ref_name
        commit = parts(1);
        return;
    end
end
end

function [bundle_hash, included_files] = hash_production_source(app_root)
files = strings(0,1);
recursive_roots = {'backend','functions','tools'};
for idx = 1:numel(recursive_roots)
    listing = dir(fullfile(app_root, recursive_roots{idx}, '**', '*.m'));
    for file_idx = 1:numel(listing)
        if ~listing(file_idx).isdir
            files(end+1,1) = string(fullfile(listing(file_idx).folder, listing(file_idx).name)); %#ok<AGROW>
        end
    end
end

root_entries = {'AAA_startup.m','launch_AAA.m','Record_analysis3.m', ...
    'AP_Analysis_APP.mlapp','VERSION','_external/DEPENDENCY_HASHES.txt'};
for idx = 1:numel(root_entries)
    candidate = fullfile(app_root, strrep(root_entries{idx}, '/', filesep));
    if isfile(candidate)
        files(end+1,1) = string(candidate); %#ok<AGROW>
    end
end

files = unique(files, 'stable');
included_files = strings(size(files));
prefix_length = strlength(string(app_root)) + 1;
for idx = 1:numel(files)
    relative = extractAfter(files(idx), prefix_length);
    included_files(idx) = replace(relative, '\', '/');
end
[included_files, order] = sort(included_files);
files = files(order);

digest = java.security.MessageDigest.getInstance('SHA-256');
for idx = 1:numel(files)
    update_digest(digest, unicode2native(char(included_files(idx)), 'UTF-8'));
    update_digest(digest, uint8(0));
    fid = fopen(files(idx), 'rb');
    if fid < 0
        error('AAA:ProvenanceReadFailed', 'Cannot read source file: %s', files(idx));
    end
    cleanup = onCleanup(@() fclose(fid));
    bytes = fread(fid, Inf, '*uint8');
    update_digest(digest, bytes);
    clear cleanup;
    update_digest(digest, uint8(10));
end
bundle_hash = digest_to_string(digest.digest());
end

function update_digest(digest, bytes)
if isempty(bytes)
    return;
end
digest.update(typecast(uint8(bytes(:)), 'int8'));
end

function value = digest_to_string(digest_bytes)
unsigned = typecast(int8(digest_bytes), 'uint8');
value = lower(string(reshape(dec2hex(unsigned, 2).', 1, [])));
end

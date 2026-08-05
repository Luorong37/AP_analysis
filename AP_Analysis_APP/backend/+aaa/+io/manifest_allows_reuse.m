function tf = manifest_allows_reuse(result_path)
%MANIFEST_ALLOWS_REUSE Reject tracked results that did not complete.
% Legacy result folders without an AAA manifest remain reusable. When an
% AAA manifest exists, its status is authoritative so a failed or pending
% run is never treated as complete merely because it left output files.

result_path = string(result_path);
if ~isscalar(result_path) || strlength(strtrim(result_path)) == 0
    tf = false;
    return;
end
if isfolder(result_path)
    manifest_file = fullfile(result_path, 'analysis_manifest.mat');
else
    manifest_file = result_path;
end
if ~isfile(manifest_file)
    tf = true;
    return;
end

try
    manifest = aaa.io.load_manifest(manifest_file);
    tf = isfield(manifest, 'status') ...
        && strcmpi(string(manifest.status), "completed");
catch
    tf = false;
end
end

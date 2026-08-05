function manifest = update_manifest(manifest_or_path, patch)
%UPDATE_MANIFEST Recursively merge a patch and rewrite the manifest.

if nargin < 2 || ~isstruct(patch) || ~isscalar(patch)
    error('AAA:IO:InvalidManifestPatch', 'patch must be one scalar struct.');
end

manifest_file = "";
if isstruct(manifest_or_path)
    manifest = manifest_or_path;
    if isfield(manifest, 'manifest_file')
        manifest_file = string(manifest.manifest_file);
    end
else
    candidate = string(manifest_or_path);
    if isfolder(candidate)
        manifest_file = string(fullfile(candidate, 'analysis_manifest.mat'));
    else
        manifest_file = candidate;
    end
    manifest = aaa.io.load_manifest(manifest_file);
end

manifest = merge_struct(manifest, patch);
if strlength(manifest_file) == 0 && isfield(manifest, 'manifest_file')
    manifest_file = string(manifest.manifest_file);
end
if strlength(manifest_file) > 0
    manifest.manifest_file = manifest_file;
    aaa.io.save_manifest(manifest, manifest_file);
end
end

function output = merge_struct(base, patch)
output = base;
names = fieldnames(patch);
for idx = 1:numel(names)
    name = names{idx};
    if isfield(output, name) && isstruct(output.(name)) && isscalar(output.(name)) ...
            && isstruct(patch.(name)) && isscalar(patch.(name))
        output.(name) = merge_struct(output.(name), patch.(name));
    else
        output.(name) = patch.(name);
    end
end
end

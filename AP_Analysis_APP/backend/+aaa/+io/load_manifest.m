function manifest = load_manifest(path_or_folder)
%LOAD_MANIFEST Load analysis_manifest.mat from a file or result folder.

path_or_folder = string(path_or_folder);
if ~isscalar(path_or_folder) || strlength(strtrim(path_or_folder)) == 0
    error('AAA:IO:InvalidManifestPath', 'Manifest path must be nonempty.');
end
if isfolder(path_or_folder)
    manifest_file = fullfile(path_or_folder, 'analysis_manifest.mat');
else
    manifest_file = path_or_folder;
end
if ~isfile(manifest_file)
    error('AAA:IO:ManifestNotFound', 'Manifest file not found: %s', manifest_file);
end
loaded = load(manifest_file, 'manifest');
if ~isfield(loaded, 'manifest') || ~isstruct(loaded.manifest) || ~isscalar(loaded.manifest)
    error('AAA:IO:InvalidManifestFile', ...
        'File does not contain one scalar manifest struct: %s', manifest_file);
end
manifest = loaded.manifest;
end

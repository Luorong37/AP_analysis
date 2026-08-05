function save_manifest(manifest, manifest_file)
%SAVE_MANIFEST Atomically replace one analysis_manifest.mat file.

if ~isstruct(manifest) || ~isscalar(manifest)
    error('AAA:IO:InvalidManifest', 'manifest must be one scalar struct.');
end
manifest_file = string(manifest_file);
if ~isscalar(manifest_file) || strlength(strtrim(manifest_file)) == 0
    error('AAA:IO:InvalidManifestPath', 'manifest_file must be nonempty.');
end
[folder,~,~] = fileparts(manifest_file);
if ~isfolder(folder)
    mkdir(folder);
end
temporary_file = string([tempname(folder), '.mat']);
cleanup = onCleanup(@() delete_if_present(temporary_file));
save(temporary_file, 'manifest', '-v7.3');
[ok, message] = movefile(temporary_file, manifest_file, 'f');
if ~ok
    error('AAA:IO:ManifestMoveFailed', ...
        'Could not replace manifest %s: %s', manifest_file, message);
end
clear cleanup;
end

function delete_if_present(path)
if isfile(path)
    delete(path);
end
end

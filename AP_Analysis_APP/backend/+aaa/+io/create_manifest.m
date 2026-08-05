function [manifest, manifest_file] = create_manifest(request, execution_plan, output_path)
%CREATE_MANIFEST Create and save analysis_manifest.mat in the output folder.

if nargin < 3 || strlength(strtrim(string(output_path))) == 0
    output_path = request.workflow.output_path;
end
output_path = strtrim(string(output_path));
if ~isscalar(output_path) || strlength(output_path) == 0
    error('AAA:IO:MissingOutputPath', ...
        'create_manifest requires one nonempty output folder.');
end
if ~isfolder(output_path)
    mkdir(output_path);
end

manifest = aaa.schema.manifest_schema(request, execution_plan);
manifest.output_path = output_path;
write_log = true;
if isfield(request, 'workflow') && isfield(request.workflow, 'write_analysis_log')
    write_log = logical(request.workflow.write_analysis_log);
end
if write_log
    manifest.log_file = string(fullfile(output_path, 'analysis.log'));
else
    manifest.log_file = "";
end
manifest_file = string(fullfile(output_path, 'analysis_manifest.mat'));
manifest.manifest_file = manifest_file;

write_manifest = true;
if isfield(request, 'workflow') && isfield(request.workflow, 'write_manifest')
    write_manifest = logical(request.workflow.write_manifest);
end
if write_manifest
    aaa.io.save_manifest(manifest, manifest_file);
else
    manifest_file = "";
    manifest.manifest_file = "";
end
end

function [ctx,outcome] = export(ctx,action)
%EXPORT Write the final mode-specific result index without copying movies.
action = lower(strtrim(string(action)));
if ~ismember(action,["run","skip"])
    error('AAA:Sections:InvalidExportAction','Export action must be run or skip.');
end
roles = channel_roles(ctx);
outcome = struct('section',"export",'action',action,'state',"completed", ...
    'channels',roles,'message',"");
if action == "skip"
    outcome.state = "skipped"; outcome.message = "Export skipped."; return;
end

% Load inputs: collect final roles, paths, request schema, and profiles.
output = string(ctx.output_path);
if strlength(output)==0, error('AAA:Sections:MissingOutputPath','Export requires output_path.'); end
if ~isfolder(output), mkdir(output); end
aaa.sections.emit_phase(ctx,"export","persist","running", ...
    "Preparing the unified final result index.");

% Run algorithm: assemble the lightweight final index without movie arrays.
index = struct('mode',ctx.mode,'scope',ctx.scope,'input_path',ctx.input_path, ...
    'output_path',output,'roles',roles,'request_schema_version', ...
    string(ctx.request.schema_version),'created_at',datetime("now"), ...
    'payload_rule',"Scientific arrays remain in channel/result section files.");
index.request=ctx.request;
index.channel_profiles=[ctx.channels.profile];
if isfield(ctx,'input_layout'),index.input_layout=ctx.input_layout;end

% Return results: leave persistence to save_analysis_results in run_plan.
ctx.shared.results.export=index;
ctx.shared.output_files.export=strings(0,1);
outcome.message = "Wrote final result index.";
end
function roles = channel_roles(ctx)
profiles = [ctx.channels.profile]; roles = string({profiles.role})';
end

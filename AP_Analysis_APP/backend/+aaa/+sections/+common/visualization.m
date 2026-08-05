function [ctx,outcome] = visualization(ctx,action)
%VISUALIZATION Generate optional overview figures from saved trace stages.
action = lower(strtrim(string(action)));
outcome = struct('section',"visualization",'action',action, ...
    'state',"completed",'channels',channel_roles(ctx),'message',"");
if action == "skip"
    outcome.state = "skipped"; outcome.message = "Visualization skipped."; return;
end

% Load inputs: reuse output references or use current in-memory results.
if action == "reuse"
    source = string(ctx.request.workflow.source_results_path);
    if strlength(source) == 0 || ~isfolder(source)
        error('AAA:Sections:VisualizationReuseSourceMissing', ...
            'visualization="reuse" requires a source result folder.');
    end
    [loaded,compatibility_report]=load_analysis_results(source,struct( ...
        'mode',ctx.mode,'scope',ctx.scope,'runtime',ctx.runtime));
    if isfield(loaded,'visualization')
        ctx.shared.output_files.visualization=loaded.visualization;
    else
        ctx.shared.output_files.visualization=struct('source_results',source);
    end
    outcome.source_files=compatibility_report.absolute_source_files;
    outcome.compatibility_report=compatibility_report;
    outcome.message = "Visualization files remain referenced from the source result.";
    return;
end
output = require_output(ctx);

% Run algorithm: render the configured analysis overview.
aaa.sections.emit_phase(ctx,"visualization","render","running", ...
    "Rendering analysis overview figures.");
spec=struct('mode',ctx.mode);
if ctx.mode=="dual"
    spec.voltage_polarity=ctx.request.params.dual.voltage_polarity;
    spec.calcium_polarity=ctx.request.params.dual.calcium_polarity;
end
output_files=render_analysis_overview(ctx.channels,output,spec,ctx.shared.results);

% Return results: register the generated figure/image bundle.
ctx.shared.output_files.visualization=output_files;
outcome.message = "Saved the analysis overview output bundle.";
end

function roles = channel_roles(ctx)
profiles = [ctx.channels.profile]; roles = string({profiles.role})';
end
function output = require_output(ctx)
output = string(ctx.output_path);
if strlength(output)==0, error('AAA:Sections:MissingOutputPath','Visualization requires output_path.'); end
if ~isfolder(output), mkdir(output); end
end

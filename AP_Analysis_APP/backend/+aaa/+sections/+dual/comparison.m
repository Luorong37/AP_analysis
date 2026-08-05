function [ctx,outcome] = comparison(ctx,action)
%COMPARISON Coordinate Dual overlap packages and persistence.
action = lower(strtrim(string(action)));
outcome = struct('section',"comparison",'action',action,'state',"completed", ...
    'channels',["voltage";"calcium"],'message',"");
if action == "skip"
    outcome.state="skipped"; outcome.message="Dual comparison skipped."; return;
end
output = require_output(ctx);

% Load inputs: reuse a comparison or collect current voltage/calcium stages.
if action == "reuse"
    aaa.sections.emit_phase(ctx,"comparison","reuse_load","running", ...
        "Loading saved Dual comparison.");
    source = string(ctx.request.workflow.source_results_path);
    [loaded,compatibility_report]=load_analysis_results(source,struct( ...
        'mode',ctx.mode,'scope',ctx.scope,'runtime',ctx.runtime));
    if ~isfield(loaded,'comparison')
        error('AAA:Sections:ComparisonReuseMissing', ...
            'Source results lack comparison.');
    end
    result = loaded.comparison;
    outcome.source_files=compatibility_report.absolute_source_files;
    outcome.compatibility_report=compatibility_report;
    message = "Reused Dual comparison.";
else
    % Run algorithm: compare current role traces and render the figures.
    aaa.sections.emit_phase(ctx,"comparison","compare","running", ...
        "Building voltage/calcium overlap packages.");
    vi = role_index(ctx,"voltage"); ci = role_index(ctx,"calcium");
    spec = struct('voltage_polarity',ctx.request.params.dual.voltage_polarity, ...
        'calcium_polarity',ctx.request.params.dual.calcium_polarity);
    result = compare_channel_traces(ctx.channels(vi).results, ...
        ctx.channels(ci).results,ctx.channels(vi).profile, ...
        ctx.channels(ci).profile,spec);
    windows=struct();if isfield(ctx.shared.results,'stim')&&isfield(ctx.shared.results.stim,'windows'),windows=ctx.shared.results.stim.windows;end
    aaa.sections.emit_phase(ctx,"comparison","render","running", ...
        "Rendering Dual comparison figures.");
    result.output_files=render_dual_comparison(result,output,windows);
    ctx.shared.output_files.comparison=result.output_files;
    message = "Built Dual overlap comparison packages.";
end

% Return results: register the comparison and output-file inventory.
ctx.shared.results.comparison = result;
aaa.sections.emit_phase(ctx,"comparison","persist","running", ...
    "Preparing the unified Dual comparison result.");
if ~isfield(ctx.shared.output_files,'comparison'),ctx.shared.output_files.comparison=strings(0,1);end
outcome.message = message;
end
function idx = role_index(ctx,role)
profiles=[ctx.channels.profile]; idx=find(string({profiles.role})==role,1);
if isempty(idx), error('AAA:Sections:MissingComparisonRole','Missing %s.',role); end
end
function output = require_output(ctx)
output=string(ctx.output_path); if strlength(output)==0, error('AAA:Sections:MissingOutputPath','Comparison requires output_path.'); end
if ~isfolder(output), mkdir(output); end
end

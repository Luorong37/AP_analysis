function [ctx,outcome] = stim(ctx,action)
%STIM Resolve and analyze visual stimulation for AP or Dual channels.
action=lower(strtrim(string(action)));
roles=string(arrayfun(@(x)x.profile.role,ctx.channels,'UniformOutput',false))';
outcome=struct('section',"stim",'action',action,'state',"completed", ...
    'channels',roles,'message',"");
if action=="skip",outcome.state="skipped";outcome.message="Stimulus analysis skipped.";return;end
output=require_output(ctx);

% Load inputs: read saved stimulus results or resolve current trial windows.
if action=="reuse"
    aaa.sections.emit_phase(ctx,"stim","reuse_load","running", ...
        "Loading saved visual-stim results.");
    source=string(ctx.request.workflow.source_results_path);
    [loaded,compatibility_report]=load_analysis_results(source,struct( ...
        'mode',ctx.mode,'scope',ctx.scope,'runtime',ctx.runtime));
    if ~isfield(loaded,'stim')
        error('AAA:Sections:StimReuseMissing', ...
            'No compatible stimulus result was found at %s.',source);
    end
    result=loaded.stim;message="Reused stimulus results.";
    outcome.source_files=compatibility_report.absolute_source_files;
    outcome.compatibility_report=compatibility_report;
else
    mode_params=ctx.request.params.(char(ctx.mode));stim_params=mode_params.stim;
    override=mode_params.stim_context_override;
    aaa.sections.emit_phase(ctx,"stim","resolve_context","running", ...
        "Resolving visual-stim metadata and trial windows.");
    [stim_context,windows]=resolve_visualstim(ctx.input_path,ctx.channels,override,stim_params);
    result=struct('supported',stim_context.supported,'status',"not_applicable", ...
        'info',stim_context,'windows',windows,'metrics',struct(), ...
        'response',struct(),'created_at',datetime("now"));
    if ~stim_context.supported
        if strcmpi(stim_context.recordmode,'visualstim')
            error('AAA:Sections:StimMetadataMissing', ...
                'Visual-stim analysis cannot continue: %s',stim_context.reason);
        end
        message="Input is not a supported visual-stim record.";
    else
        % Run algorithm: calculate responses/tuning, then render figures.
        result.status="completed";
        aaa.sections.emit_phase(ctx,"stim","metrics","running", ...
            "Computing trial responses and tuning metrics.");
        result=compute_metrics(result,ctx,stim_params);
        plot_spec=struct('mode',ctx.mode,'stim_shading_alpha',shading_alpha(ctx));
        aaa.sections.emit_phase(ctx,"stim","render","running", ...
            "Rendering visual-stim response figures.");
        result.output_files=render_visualstim_results(result,ctx.channels,output,plot_spec);
        ctx.shared.output_files.stim=result.output_files;
        message="Computed shared visual-stim windows, responses, and figures; voltage response="+result.voltage_response.effective;
        if strlength(result.voltage_response.fallback_reason)>0,message=message+" (fallback: "+result.voltage_response.fallback_reason+")";end
    end
end

% Return results: register stimulus data and generated output files.
aaa.sections.emit_phase(ctx,"stim","persist","running", ...
    "Preparing the unified visual-stim result.");
if ~isfield(ctx.shared.output_files,'stim'),ctx.shared.output_files.stim=strings(0,1);end
ctx.shared.results.stim=result;outcome.message=message;
if isfield(result,'voltage_response'),outcome.voltage_response=result.voltage_response;end
end

function result=compute_metrics(result,ctx,stim_params)
% Handler-local coordinator; external metric/tuning atoms receive resolved
% arrays, windows, polarity, and bound profiles rather than ctx/request.
windows=result.windows;metrics=struct('sensitivity',struct(),'snr',struct());
for idx=1:numel(ctx.channels)
    channel=ctx.channels(idx);role=char(lower(string(channel.profile.role)));polarity=role_polarity(ctx,channel,role);
    [sensitivity,stage]=optional_stage(channel.results,sensitivity_stages(role));
    if ~isempty(sensitivity)
        metrics.sensitivity.(role)=analyze_stimulus_trials(sensitivity,channel.profile,windows.(role),polarity);
        metrics.sensitivity.(role).input_stage=stage;
    end
    [snr,stage]=optional_stage(channel.results,snr_stages(role));
    if ~isempty(snr)
        metrics.snr.(role)=analyze_stimulus_trials(snr,channel.profile,windows.(role),polarity);
        metrics.snr.(role).input_stage=stage;
    end
end
if isempty(fieldnames(metrics.sensitivity)),metrics=rmfield(metrics,'sensitivity');end
if isempty(fieldnames(metrics.snr)),metrics=rmfield(metrics,'snr');end
result.metrics=metrics;

vi=find_role(ctx,'voltage');channel=ctx.channels(vi);requested=lower(string(stim_params.voltage_response_method));
[has_peaks,indices]=accepted_peaks(channel.results);effective=requested;fallback="";
if requested=="auto",if has_peaks,effective="peak_count";else,effective="mean_sensitivity";end,end
if effective=="peak_count"&&~has_peaks,effective="mean_sensitivity";fallback="accepted peaks unavailable";end
if effective=="peak_count"
    response=count_stimulus_peaks(indices,windows.voltage);units=response.units;
else
    if ~isfield(metrics,'sensitivity')||~isfield(metrics.sensitivity,'voltage')
        error('AAA:Sections:StimSensitivityMissing','Voltage mean-sensitivity response requires a sensitivity trace stage.');
    end
    m=metrics.sensitivity.voltage;response=struct('baseline',m.baseline_mean,'stim',m.stim_mean,'delta',m.delta_mean,'units',"mean sensitivity");units=response.units;
end
result.response.voltage=response;
result.voltage_response=struct('requested',requested,'effective',effective, ...
    'fallback_reason',fallback,'units',units,'accepted_peaks_available',has_peaks);
if windows.is_grating
    result.analysis_kind="grating_tuning";result.tuning=struct();
    result.tuning.voltage=analyze_grating_tuning(response.stim,response.baseline,windows.orientations);
    if isfield(metrics,'sensitivity')&&isfield(metrics.sensitivity,'calcium')
        m=metrics.sensitivity.calcium;result.response.calcium=struct('baseline',m.baseline_mean,'stim',m.stim_mean,'delta',m.delta_mean,'units',"mean sensitivity");
        result.tuning.calcium=analyze_grating_tuning(m.stim_mean,m.baseline_mean,windows.orientations);
    end
else
    result.analysis_kind="condition_response";
end
end

function value=role_polarity(ctx,channel,role)
if ctx.mode=="dual",value=ctx.request.params.dual.(role+"_polarity");return;end
[has_peaks,~]=accepted_peaks(channel.results);
if has_peaks,value=cell2mat(channel.results.peak_results.accepted_for_events.data.polarity(:)');return;end
[trace,~]=resolve_trace_stage(channel.results,{'sensitivity'});value=ones(1,size(trace,2));mode=lower(string(channel.profile.peak.polarity_mode));
if mode=="negative",value(:)=-1;elseif mode=="auto",for roi=1:size(trace,2),lo=abs(prctile(trace(:,roi),.5));hi=abs(prctile(trace(:,roi),99.5));if lo>hi,value(roi)=-1;end,end,end
end
function [tf,indices]=accepted_peaks(results)
tf=false;indices={};if ~isfield(results,'peak_results'),return;end;p=results.peak_results;if ~isfield(p,'accepted_for_events')||~isfield(p.accepted_for_events,'data')||~isfield(p.accepted_for_events.data,'index'),return;end
indices=p.accepted_for_events.data.index;tf=iscell(indices);
end
function [trace,stage]=optional_stage(results,names)
trace=[];stage="";try,[trace,stage]=resolve_trace_stage(results,names);catch,end
end
function names=sensitivity_stages(role),if strcmp(role,'calcium'),names={'sensitivity_smoothed','sensitivity'};else,names={'sensitivity'};end,end
function names=snr_stages(role),if strcmp(role,'calcium'),names={'snr_smoothed','snr'};else,names={'snr'};end,end
function idx=find_role(ctx,role),idx=find(arrayfun(@(x)strcmpi(string(x.profile.role),role),ctx.channels),1);end
function value=shading_alpha(ctx),if ctx.mode=="dual",value=ctx.request.params.dual.visualization.stim_shading_alpha;else,value=ctx.request.params.ap.visualization.stim_shading_alpha;end,end
function output=require_output(ctx),output=string(ctx.output_path);if strlength(output)==0,error('AAA:Sections:MissingOutputPath','Stim requires output_path.');end;if ~isfolder(output),mkdir(output);end,end

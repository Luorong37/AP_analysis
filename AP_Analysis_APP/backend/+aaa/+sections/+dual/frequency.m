function [ctx,outcome] = frequency(ctx,action)
%FREQUENCY Run/reuse Dual direct FFT and representative CWT summaries.
action=lower(strtrim(string(action)));
outcome=struct('section',"frequency",'action',action,'state',"completed", ...
    'channels',["voltage";"calcium"],'message',"");
if action=="skip", outcome.state="skipped"; outcome.message="Frequency skipped."; return; end
output=require_output(ctx);

% Load inputs: reuse frequency results or resolve current trace stages/spec.
if action=="reuse"
    aaa.sections.emit_phase(ctx,"frequency","reuse_load","running", ...
        "Loading saved time-frequency results.");
    source=string(ctx.request.workflow.source_results_path);
    [loaded,compatibility_report]=load_analysis_results(source,struct( ...
        'mode',ctx.mode,'scope',ctx.scope,'runtime',ctx.runtime));
    if ~isfield(loaded,'frequency'),error('AAA:Sections:FrequencyReuseMissing', ...
            'Source results lack frequency data.');end
    results=loaded.frequency;
    outcome.source_files=compatibility_report.absolute_source_files;
    outcome.compatibility_report=compatibility_report;
    message="Reused time-frequency results.";
else
    params=ctx.request.params.dual.frequency;
    spec=struct('min_freq_hz',params.min_freq_hz, ...
        'max_freq_hz',params.max_freq_hz, ...
        'wavelet_voices_per_octave',params.wavelet_voices_per_octave, ...
        'wavelet_name',string(params.wavelet_name));
    results=struct('parameters',spec);
    % Run algorithm: analyze each role and render the frequency figures.
    for idx=1:numel(ctx.channels)
        role=char(lower(string(ctx.channels(idx).profile.role)));
        aaa.sections.emit_phase(ctx,"frequency","fft_cwt","running", ...
            "Computing FFT and representative CWT.",string(role));
        [trace,stage]=resolve_trace_stage(ctx.channels(idx).results, ...
            {'snr','sensitivity','bleach_removed'});
        trace=double(ctx.request.params.dual.(char(role+"_polarity")))*double(trace);
        results.(role)=analyze_frequency(trace,ctx.channels(idx).profile,spec);
        results.(role).input_stage=string(stage);
    end
    windows=struct();if isfield(ctx.shared.results,'stim')&&isfield(ctx.shared.results.stim,'windows'),windows=ctx.shared.results.stim.windows;end
    aaa.sections.emit_phase(ctx,"frequency","render","running", ...
        "Rendering time-frequency figures.");
    results.output_files=render_frequency_results(results,output,windows);
    ctx.shared.output_files.frequency=results.output_files;
    message="Computed direct FFT and representative CWT summaries.";
end

% Return results: register frequency data and generated output files.
aaa.sections.emit_phase(ctx,"frequency","persist","running", ...
    "Preparing the unified time-frequency result.");
if ~isfield(ctx.shared.output_files,'frequency'),ctx.shared.output_files.frequency=strings(0,1);end
ctx.shared.results.frequency=results; outcome.message=message;
end
function output=require_output(ctx)
output=string(ctx.output_path); if strlength(output)==0, error('AAA:Sections:MissingOutputPath','Frequency requires output_path.'); end
if ~isfolder(output),mkdir(output);end
end

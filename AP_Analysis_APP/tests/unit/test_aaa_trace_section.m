function tests = test_aaa_trace_section
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
unit_dir = fileparts(mfilename('fullpath'));
app_dir = fileparts(fileparts(unit_dir));
old_path = path;
restoredefaultpath;
addpath(app_dir);
AAA_startup();
testCase.TestData.oldPath = old_path;
testCase.TestData.appDir = app_dir;
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testDualTraceRunPreservesStagesAndFiles(testCase)
[ctx, output_path, cleanup] = make_context(testCase); %#ok<ASGLU>
[ctx, outcome] = aaa.sections.common.trace(ctx, "run");
[ctx,~] = save_analysis_results(ctx,"trace",outcome);
verifyEqual(testCase, outcome.state, "completed");
verifyTrue(testCase, isfile(fullfile(output_path,'results','roi.mat')));
verifyTrue(testCase, isfile(fullfile(output_path,'results','trace.mat')));
verifyFalse(testCase, isfile(fullfile(output_path,'voltage_results.mat')));

profiles = [ctx.channels.profile];
voltage_idx = find(string({profiles.role}) == "voltage", 1);
calcium_idx = find(string({profiles.role}) == "calcium", 1);
verifyEqual(testCase, ...
    ctx.channels(voltage_idx).results.trace_results.snr.info.parent_results, ...
    {'bleach_removed','noise'});
verifyEqual(testCase, ...
    ctx.channels(calcium_idx).results.trace_results.raw_smoothed.info.parent_results, ...
    {'raw'});
end

function testTraceReuseLoadsCompleteStagesWithoutRecompute(testCase)
[source_ctx, source_path, source_cleanup] = make_context(testCase); %#ok<ASGLU>
[source_ctx, ~] = aaa.sections.common.trace(source_ctx, "run");
[source_ctx,~] = save_analysis_results(source_ctx,"trace", ...
    struct('state',"completed",'action',"run"));
[ctx, output_path, cleanup] = make_context(testCase); %#ok<ASGLU>
ctx.request.workflow.source_results_path = string(source_path);
[ctx, outcome] = aaa.sections.common.trace(ctx, "reuse");
verifyTrue(testCase,contains(outcome.message, ...
    "Saved frame rates are authoritative"));
for idx = 1:numel(ctx.channels)
    verifyTrue(testCase,isfield(ctx.channels(idx).output_files,'trace'));
    verifyEqual(testCase, ...
        fetch_trace_stage(ctx.channels(idx).results, 'sensitivity'), ...
        fetch_trace_stage(source_ctx.channels(idx).results, 'sensitivity'));
end
[ctx,~]=save_analysis_results(ctx,"trace",outcome);
verifyTrue(testCase,isfile(fullfile(output_path,'results','trace.mat')));
end

function [ctx, output_path, cleanup] = make_context(testCase)
output_path = tempname(testCase.TestData.appDir);
mkdir(output_path);
cleanup = onCleanup(@() remove_folder(output_path));
request = aaa.schema.default_request("dual", "single");
request.input_path = testCase.TestData.appDir;
request.workflow.output_path = string(output_path);
request.params.common.run_background_removal = false;
request.params.dual.bleach_mode_voltage = "linear";
request.params.dual.bleach_mode_calcium = "linear";
request.params.dual.calcium_smoothing_window = 20;
ctx = aaa.sections.create_context(request, struct());
nframes = 2048;
t = (1:nframes)';
for idx = 1:numel(ctx.channels)
    role = string(ctx.channels(idx).profile.role);
    if role == "voltage"
        raw = [100 + 0.01*t + sin(t/7), 110 + 0.02*t + cos(t/11)];
    else
        raw = [200 + 0.02*t + sin(t/31), 210 + 0.01*t + cos(t/37)];
    end
    movie_info = struct( ...
        'source_path', "synthetic", ...
        'frame_rate', ctx.channels(idx).profile.frame_rate, ...
        'frame_count', nframes, ...
        'analysis_frame_size', [4 4], ...
        'transpose_before_analysis', false, ...
        'motion', struct('applied', false));
    results = struct('movie_info',movie_info,'trace_results',struct());
    results = store_trace_stage(results, 'raw', raw, {}, ...
        "roi.mat", movie_info, 'synthetic_roi_mean', struct());
    ctx.channels(idx).results = results;
end
end

function remove_folder(path_value)
if isfolder(path_value)
    rmdir(path_value, 's');
end
end

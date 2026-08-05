function tests = test_aaa_peak_section
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

function testDualRunPersistsModernPeakChain(testCase)
[ctx, output_path, cleanup] = make_context(testCase, "dual"); %#ok<ASGLU>
[ctx, outcome] = aaa.sections.common.peak(ctx, "run");
[ctx,~]=save_analysis_results(ctx,"peak",outcome);
verifyEqual(testCase, outcome.state, "completed");
verifyTrue(testCase,isfile(fullfile(output_path,'results','peak.mat')));
verifyTrue(testCase, isfield(ctx.channels(1).results.peak_results, ...
    'accepted_for_events'));
verifyEqual(testCase, ...
    ctx.channels(1).results.peak_results.current_result, ...
    'accepted_for_events');
end

function testReuseIsStrictAndUsesCurrentTraceShape(testCase)
[source_ctx, source_path, source_cleanup] = ...
    make_context(testCase, "dual"); %#ok<ASGLU>
[source_ctx, ~] = aaa.sections.common.peak(source_ctx, "run");
[source_ctx,~]=save_analysis_results(source_ctx,"peak", ...
    struct('state',"completed",'action',"run"));

[ctx, output_path, cleanup] = make_context(testCase, "dual"); %#ok<ASGLU>
ctx.request.workflow.source_results_path = string(source_path);
[ctx, outcome] = aaa.sections.common.peak(ctx, "reuse");
[ctx,~]=save_analysis_results(ctx,"peak",outcome);
verifyEqual(testCase, outcome.state, "completed");
verifyTrue(testCase, ...
    ctx.channels(1).results.peak_results.reuse_info.reused);
verifyTrue(testCase,isfile(fullfile(output_path,'results','peak.mat')));

ctx.channels(1).results.trace_results.sensitivity.data = zeros(5, 1);
verifyError(testCase, ...
    @() aaa.sections.common.peak(ctx, "reuse"), ...
    'AAA:Sections:InvalidPeakReuse');
end

function [ctx, output_path, cleanup] = make_context(testCase, mode)
output_path = tempname(testCase.TestData.appDir);
mkdir(output_path);
cleanup = onCleanup(@() remove_folder(output_path));
request = aaa.schema.default_request(mode, "single");
request.input_path = testCase.TestData.appDir;
request.workflow.output_path = string(output_path);
if mode == "dual"
    request.params.dual.run_manual_voltage_peak_edit = false;
else
    request.params.ap.run_manual_peak_refinement = false;
end
ctx = aaa.sections.create_context(request, struct());
profiles = [ctx.channels.profile];
voltage_idx = find(string({profiles.role}) == "voltage", 1);
t = (1:120)';
sensitivity = sin(2*pi*t/20);
movie_info = struct( ...
    'source_path', "synthetic", ...
    'transpose_before_analysis', false, ...
    'motion', struct('applied', false));
ctx.channels(voltage_idx).results = store_trace_stage( ...
    struct(), 'sensitivity', sensitivity, ...
    {'bleach_removed','baseline'}, "roi.mat", movie_info, ...
    'ratio', struct());
end

function remove_folder(path_value)
if isfolder(path_value)
    rmdir(path_value, 's');
end
end

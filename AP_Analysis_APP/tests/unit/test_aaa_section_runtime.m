function tests = test_aaa_section_runtime
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
unit_dir = fileparts(mfilename('fullpath'));
app_dir = fileparts(fileparts(unit_dir));
backend_dir = fullfile(app_dir, 'backend');
addpath(backend_dir);
testCase.TestData.backendDir = backend_dir;
testCase.TestData.tempRoot = tempname;
mkdir(testCase.TestData.tempRoot);
end

function teardownOnce(testCase)
rmpath(testCase.TestData.backendDir);
if isfolder(testCase.TestData.tempRoot)
    rmdir(testCase.TestData.tempRoot, 's');
end
end

function testContextUsesCanonicalChannelProfiles(testCase)
ap_request = aaa.schema.default_request("ap", "single");
ap_request.input_path = string(testCase.TestData.tempRoot);
ap_ctx = aaa.sections.create_context(ap_request, struct());
verifyNumElements(testCase, ap_ctx.channels, 1);
verifyEqual(testCase, ap_ctx.channels(1).profile.role, "voltage");
verifyEqual(testCase, ap_ctx.channels(1).profile.frame_rate, 400);
verifyFalse(testCase, isfield(ap_ctx.channels, 'role'));
verifyFalse(testCase, isfield(ap_ctx.channels, 'frame_rate'));
verifyEqual(testCase, ...
    ap_ctx.channels(1).profile.trace.noise_reference.method, "wdenoise");

dual_request = aaa.schema.default_request("dual", "single");
dual_request.input_path = string(testCase.TestData.tempRoot);
dual_ctx = aaa.sections.create_context(dual_request, struct());
verifyNumElements(testCase, dual_ctx.channels, 2);
profiles = [dual_ctx.channels.profile];
verifyEqual(testCase, string({profiles.role}), ["voltage", "calcium"]);
verifyEqual(testCase, [dual_ctx.channels.camera_index], [2 1]);
verifyEqual(testCase, ...
    dual_ctx.channels(2).profile.trace.noise_reference.normalized_cutoff, ...
    0.05, 'AbsTol', eps);
verifyEqual(testCase, ...
    dual_ctx.channels(2).profile.trace.smoothing.method, "movmean");
end

function testRunReuseSkipEmitExactEvents(testCase)
received = cell(0, 1);
runtime = struct('notify', @capture_event);
request = aaa.schema.default_request("dual", "single");
request.input_path = string(testCase.TestData.tempRoot);
request.workflow.output_path = string(fullfile(testCase.TestData.tempRoot,'runtime_events'));
ctx = aaa.sections.create_context(request, runtime);
overrides = struct( ...
    'input', @successful_handler, ...
    'motion', @successful_handler, ...
    'registration', @successful_handler);
entries = aaa.sections.registry("dual", overrides);
plan = struct('input',"run",'motion',"reuse",'registration',"skip");

[ctx, receipt] = aaa.sections.run_plan(ctx, plan, entries);

verifyEqual(testCase, receipt.status, "completed");
verifyEqual(testCase, receipt.completed_sections, "input");
verifyEqual(testCase, receipt.reused_sections, "motion");
verifyEqual(testCase, receipt.skipped_sections, "registration");
events = vertcat(received{:});
events = events(string({events.type})=="section");
verifyEqual(testCase, string({events.state}), ...
    ["running","completed","reusing","reused","skipped"]);
verifyEqual(testCase, string(arrayfun( ...
    @(event) event.details.action, events, 'UniformOutput', false)), ...
    ["run";"run";"reuse";"reuse";"skip"]);
verifyEqual(testCase, ctx.channels(1).status.section, "motion");
verifyEqual(testCase, ctx.channels(1).status.state, "reused");

    function capture_event(event)
        received{end + 1, 1} = event;
    end
end

function testReportedChannelFailureStopsPlanAndReturnsReceipt(testCase)
request = aaa.schema.default_request("dual", "single");
request.input_path = string(testCase.TestData.tempRoot);
request.workflow.output_path = string(fullfile(testCase.TestData.tempRoot,'runtime_failure'));
ctx = aaa.sections.create_context(request, struct());
overrides = struct( ...
    'input', @successful_handler, ...
    'motion', @reported_failure_handler, ...
    'registration', @marks_later_execution);
entries = aaa.sections.registry("dual", overrides);
plan = struct('input',"run",'motion',"run",'registration',"run");

[ctx, receipt] = aaa.sections.run_plan(ctx, plan, entries);

verifyEqual(testCase, receipt.status, "failed");
verifyEqual(testCase, receipt.failed_section, "motion");
verifyEqual(testCase, receipt.failed_channel, "calcium");
verifyEqual(testCase, receipt.error.identifier, "AAA:Test:CalciumFailure");
verifyTrue(testCase, isa(ctx.execution.exception, 'MException'));
verifyFalse(testCase, isfield(ctx.shared.results, 'later_executed'));
verifyEqual(testCase, ctx.channels(2).status.state, "failed");
verifyEqual(testCase, receipt.section_status.motion.state, "failed");
end

function testThrownFailureIsCapturedForUpperWorkflowRethrow(testCase)
request = aaa.schema.default_request("ap", "single");
request.input_path = string(testCase.TestData.tempRoot);
request.workflow.output_path = string(fullfile(testCase.TestData.tempRoot,'runtime_throw'));
ctx = aaa.sections.create_context(request, struct());
entries = aaa.sections.registry("ap", struct('input',@throwing_handler));

[ctx, receipt] = aaa.sections.run_plan(ctx, struct('input',"run"), entries);

verifyEqual(testCase, receipt.status, "failed");
verifyEqual(testCase, receipt.failed_section, "input");
verifyEqual(testCase, receipt.failed_channel, "");
verifyEqual(testCase, receipt.error.identifier, "AAA:Test:ThrownFailure");
verifyEqual(testCase, ctx.execution.exception.identifier, ...
    'AAA:Test:ThrownFailure');
end

function testReceiptNeverCopiesScientificPayloads(testCase)
request = aaa.schema.default_request("ap", "single");
request.input_path = string(testCase.TestData.tempRoot);
ctx = aaa.sections.create_context(request, struct());
ctx.channels(1).data.movie_3d = rand(12, 9, 5);
ctx.channels(1).results.trace = rand(100, 4);
ctx.channels(1).output_files.debug_payload = rand(20);

receipt = aaa.sections.build_receipt(ctx);

verifyFalse(testCase, isfield(receipt, 'request'));
verifyFalse(testCase, isfield(receipt, 'channels'));
verifyFalse(testCase, isfield(receipt.channel_status, 'data'));
verifyFalse(testCase, isfield(receipt.channel_status, 'results'));
verifyFalse(testCase, isfield(receipt.channel_status, 'output_files'));
verifyEqual(testCase, receipt.channel_status.role, "voltage");
end

function [ctx, outcome] = successful_handler(ctx, action)
roles = channel_roles(ctx);
ctx.shared.results.last_successful_action = string(action);
outcome = struct( ...
    'state', "completed", ...
    'channels', roles, ...
    'message', "Test handler completed.");
end

function [ctx, outcome] = reported_failure_handler(ctx, ~)
outcome = struct( ...
    'state', "failed", ...
    'channels', ["voltage";"calcium"], ...
    'failed_channel', "calcium", ...
    'message', "Synthetic calcium failure.", ...
    'error', struct( ...
        'identifier',"AAA:Test:CalciumFailure", ...
        'message',"Synthetic calcium failure."));
end

function [ctx, outcome] = marks_later_execution(ctx, ~)
ctx.shared.results.later_executed = true;
outcome = struct('state',"completed");
end

function [ctx, outcome] = throwing_handler(ctx, ~) %#ok<INUSD>
outcome = struct(); %#ok<NASGU>
error('AAA:Test:ThrownFailure', 'Synthetic thrown failure.');
end

function roles = channel_roles(ctx)
profiles = [ctx.channels.profile];
roles = string({profiles.role})';
end

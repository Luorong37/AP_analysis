function tests = test_aaa_io
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

function testManifestRoundTrip(testCase)
request = aaa.schema.default_request("ap", "single");
request.input_path = string(testCase.TestData.tempRoot);
output_path = fullfile(testCase.TestData.tempRoot, 'result');
request.workflow.output_path = string(output_path);
[~,~,~,plan] = aaa.schema.resolve_sections("ap", "full", struct());
[manifest, manifest_file] = aaa.io.create_manifest(request, plan, output_path);
verifyTrue(testCase, isfile(manifest_file));
verifyEqual(testCase, manifest.status, "pending");
verifyEqual(testCase, manifest.source.app_version, "0.1.0");
verifyEqual(testCase, manifest.source.hashes.algorithm, "SHA-256");
verifyGreaterThan(testCase, strlength(manifest.source.hashes.source_bundle), 0);
verifyGreaterThan(testCase, manifest.source.hashes.file_count, 0);
verifyGreaterThan(testCase, strlength(manifest.source.matlab_release), 0);

patch = struct('status',"completed",'completed_sections',string(fieldnames(plan)), ...
    'completed_at',datetime('now'),'output_files',"trace_results.mat");
aaa.io.update_manifest(manifest_file, patch);
loaded = aaa.io.load_manifest(output_path);
verifyEqual(testCase, loaded.status, "completed");
verifyEqual(testCase, loaded.output_files, "trace_results.mat");
verifyEqual(testCase, loaded.manifest_file, string(manifest_file));
end

function testLoggerAndEventCallback(testCase)
log_file = fullfile(testCase.TestData.tempRoot, 'event.log');
logger = aaa.io.AnalysisLogger(log_file, false);
received = struct();
runtime = struct('logger',logger,'notify',@capture_event);
event = aaa.io.emit_event(runtime, struct( ...
    'mode',"dual",'scope',"single",'section',"trace", ...
    'state',"running",'message',"Trace started"));
logger.close();
verifyEqual(testCase, event.type, "status");
verifyEqual(testCase, received.section, "trace");
verifyTrue(testCase, isfile(log_file));
verifyTrue(testCase, contains(string(fileread(log_file)), "Trace started"));

    function capture_event(value)
        received = value;
    end
end

function testManifestHonorsDisabledLog(testCase)
request = aaa.schema.default_request("dual", "single");
request.input_path = string(testCase.TestData.tempRoot);
request.workflow.write_analysis_log = false;
output_path = fullfile(testCase.TestData.tempRoot, 'no_log_result');
[~,~,~,plan] = aaa.schema.resolve_sections("dual", "full", struct());
[manifest, manifest_file] = aaa.io.create_manifest(request, plan, output_path);
verifyTrue(testCase, isfile(manifest_file));
verifyEqual(testCase, manifest.log_file, "");
end

function testAttachLoggerFansOutToResultFolder(testCase)
request = aaa.schema.default_request("ap", "single");
external_file = fullfile(testCase.TestData.tempRoot, 'external.log');
result_dir = fullfile(testCase.TestData.tempRoot, 'fanout_result');
external_logger = aaa.io.AnalysisLogger(external_file, false);
runtime = struct('logger', external_logger);
[runtime, result_logger, result_file, owns_logger] = ...
    aaa.io.attach_analysis_logger(runtime, request, result_dir);
verifyTrue(testCase, owns_logger);
verifyEqual(testCase, result_file, string(fullfile(result_dir, 'analysis.log')));
verifyFalse(testCase, result_logger.EchoToCommandWindow);
aaa.io.emit_event(runtime, struct('message',"fanout marker",'state',"running"));
result_logger.close();
external_logger.close();
verifyTrue(testCase, contains(string(fileread(result_file)), "fanout marker"));
verifyTrue(testCase, contains(string(fileread(external_file)), "fanout marker"));
end

function testManifestControlsReuseEligibility(testCase)
result_dir = fullfile(testCase.TestData.tempRoot, 'reuse_status');
request = aaa.schema.default_request("ap", "single");
request.input_path = string(testCase.TestData.tempRoot);
[~,~,~,plan] = aaa.schema.resolve_sections("ap", "full", struct());
[~, manifest_file] = aaa.io.create_manifest(request, plan, result_dir);

verifyFalse(testCase, aaa.io.manifest_allows_reuse(result_dir));
aaa.io.update_manifest(manifest_file, struct('status', "failed"));
verifyFalse(testCase, aaa.io.manifest_allows_reuse(result_dir));
aaa.io.update_manifest(manifest_file, struct('status', "completed"));
verifyTrue(testCase, aaa.io.manifest_allows_reuse(result_dir));
end

function testLegacyResultWithoutManifestRemainsReusable(testCase)
legacy_dir = fullfile(testCase.TestData.tempRoot, 'legacy_result');
mkdir(legacy_dir);
verifyTrue(testCase, aaa.io.manifest_allows_reuse(legacy_dir));
end

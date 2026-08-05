function tests = test_aaa_schema
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

function testDefaultRequestContract(testCase)
request = aaa.schema.default_request("dual", "record");
verifyEqual(testCase, request.schema_version, "1.0.0");
verifyEqual(testCase, request.mode, "dual");
verifyEqual(testCase, request.scope, "record");
verifyTrue(testCase, all(isfield(request, ...
    {'input_path','workflow','params','record'})));
verifyTrue(testCase, all(isfield(request.workflow, ...
    {'preset','sections','source_results_path','auto_detect_reuse_source', ...
    'roi_path','output_path','run_name'})));
verifyTrue(testCase, request.workflow.auto_detect_reuse_source);
verifyTrue(testCase, all(isfield(request.params, {'common','ap','dual'})));
verifyEqual(testCase, request.params.dual.voltage_frame_rate, 400);
verifyEqual(testCase, request.params.dual.calcium_transpose_movie, true);
verifyEqual(testCase, request.record.skip_existing, false);
verifyEqual(testCase, request.record.output_root, "Dual_analysis3");
end

function testEveryBooleanIsPublicCheckbox(testCase)
catalog = aaa.schema.parameter_catalog();
logical_entries = catalog([catalog.value_type] == "logical");
verifyNotEmpty(testCase, logical_entries);
verifyTrue(testCase, all([logical_entries.public]));
verifyTrue(testCase, all([logical_entries.ui_control] == "checkbox"));
verifyTrue(testCase, all(strlength([catalog.ui_level]) > 0));
verifyTrue(testCase, all(ismember([catalog.ui_level], ["Basic","Advanced","Expert"])));
end


function testBackgroundDefaultsArePublicAndShared(testCase)
ap = aaa.schema.default_request("ap", "single");
dual = aaa.schema.default_request("dual", "single");
verifyTrue(testCase, ap.params.common.run_background_removal);
verifyEqual(testCase, ap.params.common.background, ...
    dual.params.common.background);
verifyEqual(testCase, ap.params.common.background.inner_distance, 2);
verifyEqual(testCase, ap.params.common.background.outer_distance, 32);
verifyEqual(testCase, ap.params.common.background.background_fraction, 0.05);
verifyEqual(testCase, ap.params.common.background.fit_span, 0.01);
verifyEqual(testCase, ap.params.common.background.fit_method, "rloess");

catalog = aaa.schema.parameter_catalog("ap", "single");
row = catalog([catalog.key] == "params.common.run_background_removal");
verifyEqual(testCase, row.ui_control, "checkbox");
verifyTrue(testCase, row.editable);
end

function testBackgroundAnnulusValidation(testCase)
request = aaa.schema.default_request("dual", "single");
request.input_path = string(testCase.TestData.tempRoot);
request.params.common.background.outer_distance = ...
    request.params.common.background.inner_distance;
[~, report] = aaa.schema.validate_request(request);
verifyFalse(testCase, report.valid);
verifyTrue(testCase, any(contains(report.errors, ...
    "outer_distance must be greater")));
end

function testAPRetiredBackendsAndManualGatingAreRemoved(testCase)
catalog = aaa.schema.parameter_catalog("ap", "single");
keys=[catalog.key];
verifyFalse(testCase,any(contains(keys,"volpy",'IgnoreCase',true)));
verifyFalse(testCase,any(contains(keys,"nonrigid",'IgnoreCase',true)));
verifyFalse(testCase,any(keys=="params.ap.run_manual_peak_gating"));
request = aaa.schema.default_request("ap", "single");
verifyTrue(testCase,request.params.ap.run_manual_peak_refinement);
verifyEmpty(testCase, request.params.ap.peak_min_distance_frames);
end

function testAPStimulusDefaultsArePublic(testCase)
request = aaa.schema.default_request("ap", "single");
verifyEqual(testCase,request.params.ap.camera_index,1);
verifyEqual(testCase,request.params.ap.stim.voltage_response_method,"auto");
verifyEqual(testCase,request.params.ap.stim.flash_baseline_window_s,2);
catalog=aaa.schema.section_catalog("ap");
verifyTrue(testCase,any([catalog.name]=="stim"));
[~,~,~,plan]=aaa.schema.resolve_sections("ap","full",struct());
verifyEqual(testCase,plan.stim,"run");
end

function testDualPresetParity(testCase)
[preset,~,~,plan] = aaa.schema.resolve_sections("dual", "analysis_only", struct());
verifyEqual(testCase, preset, "analysis_only");
verifyEqual(testCase, plan.input, "skip");
verifyEqual(testCase, plan.trace, "reuse");
verifyEqual(testCase, plan.peak, "reuse");
verifyEqual(testCase, plan.frequency, "run");
verifyEqual(testCase, plan.export, "run");
end

function testAPStatisticsIsOneSection(testCase)
catalog = aaa.schema.section_catalog("ap");
names = [catalog.name];
verifyTrue(testCase, ismember("ap_statistics", names));
verifyFalse(testCase, any(ismember(["trend","average_ap","sequence","isi"], names)));
[~,~,~,plan] = aaa.schema.resolve_sections("ap", "analysis_only", struct());
verifyEqual(testCase, plan.ap_events, "skip");
verifyEqual(testCase, plan.ap_statistics, "skip");
end

function testCustomRequiresCompletePlan(testCase)
verifyError(testCase, ...
    @() aaa.schema.resolve_sections("ap", "custom", struct('input',"run")), ...
    'AAA:Schema:IncompleteCustomSections');
end

function testDependencyError(testCase)
overlay = struct('trace',"skip",'peak',"run");
verifyError(testCase, ...
    @() aaa.schema.resolve_sections("dual", "full", overlay), ...
    'AAA:Schema:InvalidSectionDependency');
end

function testValidationAppliesDefaults(testCase)
request = struct('mode',"ap",'scope',"single", ...
    'input_path',string(testCase.TestData.tempRoot));
[request, report] = aaa.schema.validate_request(request);
verifyTrue(testCase, report.valid, strjoin(report.errors, newline));
verifyEqual(testCase, request.workflow.preset, "full");
verifyEqual(testCase, report.execution_plan.ap_statistics, "skip");
verifyEqual(testCase, request.params.ap.AP_window_width, 15);
end

function testValidationRejectsUnknownField(testCase)
request = aaa.schema.default_request("dual", "single");
request.input_path = string(testCase.TestData.tempRoot);
request.params.dual.typo_parameter = 1;
[~, report] = aaa.schema.validate_request(request);
verifyFalse(testCase, report.valid);
verifyTrue(testCase, any(contains(report.errors, "typo_parameter")));
end

function testReuseRequiresSource(testCase)
request = aaa.schema.default_request("dual", "single");
request.input_path = string(testCase.TestData.tempRoot);
request.workflow.preset = "analysis_only";
[~, report] = aaa.schema.validate_request(request);
verifyFalse(testCase, report.valid);
verifyTrue(testCase, any(contains(report.errors, "reuse source")));
end

function testRecordReuseDefersSourceResolution(testCase)
request = aaa.schema.default_request("dual", "record");
request.input_path = string(testCase.TestData.tempRoot);
request.workflow.preset = "analysis_only";
[~, report] = aaa.schema.validate_request(request);
verifyTrue(testCase, report.valid, strjoin(report.errors, newline));
verifyTrue(testCase, any(contains(report.warnings, "scheduler")));
end

function testAPRecordValidationDoesNotRequireDualRunMotion(testCase)
request = aaa.schema.default_request("ap", "record");
request.input_path = string(testCase.TestData.tempRoot);
verifyFalse(testCase, isfield(request.record, 'run_motion'));
[~, report] = aaa.schema.validate_request(request);
verifyTrue(testCase, report.valid, strjoin(report.errors, newline));
end

function testDirectRunSpecificRequestWins(testCase)
generic = aaa.schema.default_request("dual", "single");
generic.input_path = "generic";
specific = aaa.schema.default_request("dual", "single");
specific.input_path = "specific";
assignin('base', 'aaa_request', generic);
assignin('base', 'dual3_request', specific);
cleanup = onCleanup(@clear_direct_run_requests); %#ok<NASGU>
[resolved, source_name] = aaa.helpers.common.direct_run_request("dual", "single");
verifyEqual(testCase, resolved.input_path, "specific");
verifyEqual(testCase, source_name, "dual3_request");
end

function clear_direct_run_requests()
evalin('base', 'clear aaa_request dual3_request');
end

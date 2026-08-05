function tests = test_aaa_reuse_resolution
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
unit_dir = fileparts(mfilename('fullpath'));
app_dir = fileparts(fileparts(unit_dir));
testCase.TestData.oldPath = path;
addpath(app_dir);
AAA_startup();
testCase.TestData.root = string(tempname);
mkdir(testCase.TestData.root);
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
if isfolder(testCase.TestData.root), rmdir(testCase.TestData.root,'s'); end
end

function testNewestCompletedCompatibleCandidateWins(testCase)
cycle = fullfile(testCase.TestData.root,'RecordA','Cycle1');
old = fullfile(cycle,'AP_analysis3','old');
new = fullfile(cycle,'AP_analysis3','new');
failed = fullfile(cycle,'AP_analysis3','failed_newest');
write_ap_result(old,"completed",datetime(2026,1,1));
write_ap_result(new,"completed",datetime(2026,2,1));
write_ap_result(failed,"failed",datetime(2026,3,1));
[~,~,~,plan] = aaa.schema.resolve_sections("ap","reuse_trace",struct());
[source,report] = aaa.helpers.common.resolve_reuse_source( ...
    "ap",cycle,plan,struct('cycle_name',"Cycle1",'label',"Cycle1"));
verifyEqual(testCase,source,string(new));
verifyEqual(testCase,report.status,"resolved_latest");
verifyTrue(testCase,any(~[report.candidates.valid]));
end

function testRecordRootResolvesExactCycle(testCase)
record = fullfile(testCase.TestData.root,'RecordB');
one = fullfile(record,'Cycle1','AP_analysis3','run');
two = fullfile(record,'Cycle2','AP_analysis3','run');
write_ap_result(one,"completed",datetime(2026,1,1));
write_ap_result(two,"completed",datetime(2026,2,1));
[~,~,~,plan] = aaa.schema.resolve_sections("ap","reuse_trace",struct());
source = aaa.helpers.common.resolve_reuse_source( ...
    "ap",record,plan,struct('cycle_name',"Cycle1",'label',"Cycle1"));
verifyEqual(testCase,source,string(one));
end

function testTraceReuseRestoresResultsWithoutMovie(testCase)
source = fullfile(testCase.TestData.root,'input_reuse');
mkdir(source);
trace_results = minimal_trace_results(7,2); %#ok<NASGU>
save(fullfile(source,'trace_results.mat'),'trace_results');
request = aaa.schema.default_request("ap","single");
request.input_path = testCase.TestData.root;
request.workflow.source_results_path = string(source);
request.workflow.preset="reuse_trace";
request.workflow.output_path=fullfile(testCase.TestData.root,'trace_reuse_output');
request.params.ap.freq=20;
ctx=aaa.sections.create_context(request,struct());
[ctx,outcome] = aaa.sections.common.trace(ctx,"reuse");
verifyFalse(testCase,isfield(ctx.channels(1).data,'movie_3d'));
verifyFalse(testCase,ctx.channels(1).data.movie_info.movie_loaded);
verifyEqual(testCase,ctx.channels(1).data.time,(1:7)'/20);
verifyEqual(testCase,ctx.channels(1).results.trace_results.raw.data,zeros(7,2));
verifyFalse(testCase,outcome.movie_loaded);
end

function testSemanticResultUsesContentsNotFilename(testCase)
folder=fullfile(testCase.TestData.root,'semantic','Cycle1','odd','run');
mkdir(folder);
mystery=struct('movie_info',struct('role',"voltage",'frame_rate',37), ...
    'trace_results',minimal_trace_results(6,2,37)); %#ok<NASGU>
save(fullfile(folder,'anything_results.mat'),'mystery');
[~,~,~,plan]=aaa.schema.resolve_sections("ap","reuse_trace",struct());
[source,report]=aaa.helpers.common.resolve_reuse_source("ap", ...
    fullfile(testCase.TestData.root,'semantic','Cycle1'),plan, ...
    struct('cycle_name',"Cycle1",'label',"Cycle1"));
verifyEqual(testCase,source,string(folder));
verifyEqual(testCase,report.chosen_inspection.saved_frame_rates.voltage,37);
verifyEqual(testCase,report.chosen_inspection.trace_bindings.variable_name,"result");
verifyEqual(testCase, ...
    report.chosen_inspection.compatibility_report.adapter_used, ...
    "adapt_legacy_results");
end

function testSavedFrameRateOverridesReuseRequest(testCase)
folder=fullfile(testCase.TestData.root,'fps_source');mkdir(folder);
mystery=struct('movie_info',struct('role',"voltage",'frame_rate',37), ...
    'trace_results',minimal_trace_results(6,2,37)); %#ok<NASGU>
save(fullfile(folder,'renamed_results.mat'),'mystery');
request=aaa.schema.default_request("ap","single");
request.input_path=testCase.TestData.root;
request.workflow.preset="reuse_trace";
request.workflow.source_results_path=string(folder);
request.params.ap.freq=999;
[request,resolution]=aaa.helpers.common.resolve_request_reuse_source(request);
verifyEqual(testCase,request.params.ap.freq,37);
verifyEqual(testCase,resolution.frame_rate_overrides.ui_value,999);
verifyEqual(testCase,resolution.frame_rate_overrides.saved_value,37);
end

function testSplitMovieDiscoveryAndManifestProvenance(testCase)
cycle=fullfile(testCase.TestData.root,'SplitCycle');camera=fullfile(cycle,'Cam1_voltage');
mkdir(camera);
part1=fullfile(camera,'movie_stack01.tif');
part2=fullfile(camera,'movie_stack02.tif');
imwrite(uint16(ones(3,4)),part1);imwrite(uint16(2*ones(3,4)),part2);
manifest=struct('actual',struct('movie_paths',{{fullfile(camera,'movie.tif')}}), ...
    'spec',struct('labels',{{'voltage'}})); %#ok<NASGU>
save(fullfile(cycle,'cycle_manifest.mat'),'manifest');
request=aaa.schema.default_request("ap","single");
request.input_path=string(cycle);
request.workflow.output_path=fullfile(testCase.TestData.root,'split_output');
request.workflow.preset="custom";
[~,~,~,base_plan]=aaa.schema.resolve_sections("ap","full",struct());
names=fieldnames(base_plan);
for idx=1:numel(names),request.workflow.sections.(names{idx})="skip";end
request.workflow.sections.input="run";
receipt=aaa.workflows.AP_analysis3(request,struct());
saved=aaa.io.load_manifest(receipt.manifest_file);
verifyEqual(testCase,saved.input_sources.part_count,2);
verifyTrue(testCase,saved.input_sources.is_split);
verifyEqual(testCase,string(saved.input_sources.physical_files), ...
    string({part1;part2}));
verifyEqual(testCase,saved.input_sources.logical_path,string(fullfile(camera,'movie.tif')));
end

function testEmptyRecordCycleBlocksMovieWorkflow(testCase)
record=fullfile(testCase.TestData.root,'EmptyRecord');
mkdir(fullfile(record,'Cycle1'));
request=aaa.schema.default_request("ap","record");
request.input_path=string(record);
[~,~,~,plan]=aaa.schema.resolve_sections("ap","full",struct());
inspection=aaa.helpers.common.inspect_input_structure(request,plan);
verifyFalse(testCase,inspection.valid);
verifyTrue(testCase,any(contains(inspection.errors,"empty")));
end

function testMissingMovieRequiresReusePresetAndSource(testCase)
record=fullfile(testCase.TestData.root,'NoMovieRecord');
mkdir(fullfile(record,'Cycle1','AP_analysis3','old'));
request=aaa.schema.default_request("ap","record");
request.input_path=string(record);
[~,~,~,plan]=aaa.schema.resolve_sections("ap","full",struct());
inspection=aaa.helpers.common.inspect_input_structure(request,plan);
verifyFalse(testCase,inspection.valid);
verifyTrue(testCase,any(contains(inspection.errors,"enter a manual Reuse source")));
end

function testManualReuseOverridesAutomaticCandidate(testCase)
record=fullfile(testCase.TestData.root,'PriorityRecord');
cycle=fullfile(record,'Cycle1');mkdir(cycle);
auto=fullfile(cycle,'AP_analysis3','auto');
manual_root=fullfile(testCase.TestData.root,'ManualResults');
manual=fullfile(manual_root,'Cycle1','AP_analysis3','manual');
write_ap_result(auto,"completed",datetime(2026,5,1));
write_ap_result(manual,"completed",datetime(2026,4,1));
request=aaa.schema.default_request("ap","record");
request.input_path=string(record);
request.workflow.preset="reuse_trace";
request.workflow.source_results_path=string(manual_root);
[~,~,~,plan]=aaa.schema.resolve_sections("ap","reuse_trace",struct());
inspection=aaa.helpers.common.inspect_input_structure(request,plan);
verifyTrue(testCase,inspection.valid,strjoin(inspection.errors,newline));
verifyEqual(testCase,inspection.cycles.selected_reuse_source,string(manual));
verifyEqual(testCase,inspection.cycles.auto_reuse_candidate,string(auto));
verifyTrue(testCase,inspection.cycles.manual_overrode_auto);
end

function testDisabledAutoDetectionRequiresManualSource(testCase)
record=fullfile(testCase.TestData.root,'AutoDisabledRecord');
auto=fullfile(record,'Cycle1','AP_analysis3','auto');
write_ap_result(auto,"completed",datetime(2026,6,1));
request=aaa.schema.default_request("ap","record");
request.input_path=string(record);
request.workflow.preset="reuse_trace";
request.workflow.auto_detect_reuse_source=false;
[~,~,~,plan]=aaa.schema.resolve_sections("ap","reuse_trace",struct());
inspection=aaa.helpers.common.inspect_input_structure(request,plan);
verifyFalse(testCase,inspection.valid);
verifyTrue(testCase,any(contains(inspection.errors,"disabled")));
end

function testAPReuseWorkflowCreatesNewOutputWithoutMovie(testCase)
source = fullfile(testCase.TestData.root,'workflow_source');
output = fullfile(testCase.TestData.root,'workflow_output');
input = fullfile(testCase.TestData.root,'workflow_input');
mkdir(input);
write_ap_result(source,"completed",datetime(2026,4,1));
rois = struct('bwmask',true(2)); %#ok<NASGU>
save(fullfile(source,'1_raw_ROI.mat'),'rois');
request = aaa.schema.default_request("ap","single");
request.input_path = string(input);
request.workflow.preset = "reuse_trace";
request.workflow.source_results_path = string(source);
request.workflow.output_path = string(output);
request.workflow.sections.peak = "skip";
request.workflow.sections.stim = "skip";
request.workflow.sections.visualization = "skip";
receipt = aaa.workflows.AP_analysis3(request,struct());
verifyEqual(testCase,receipt.status,"completed");
verifyTrue(testCase,isfile(fullfile(output,'results','trace.mat')));
verifyTrue(testCase,isfile(fullfile(output,'results','roi.mat')));
verifyTrue(testCase,isfile(fullfile(output,'results','export.mat')));
verifyFalse(testCase,any(isfile(fullfile(output,["movie.mat","movie.tif"]))));
end

function write_ap_result(folder,status,completed_at)
if ~isfolder(folder), mkdir(folder); end
trace_results = minimal_trace_results(5,1); %#ok<NASGU>
save(fullfile(folder,'trace_results.mat'),'trace_results');
manifest = struct('status',string(status),'completed_at',completed_at); %#ok<NASGU>
save(fullfile(folder,'analysis_manifest.mat'),'manifest');
end

function trace_results = minimal_trace_results(nframes,nrois,frame_rate)
if nargin<3,frame_rate=20;end
trace_results = struct();
for name = ["raw","bleach_removed","baseline","noise_reference", ...
        "noise","sensitivity","snr"]
    trace_results.(name) = struct('data',zeros(nframes,nrois), ...
        'frame_rate',frame_rate,'info',struct());
end
end

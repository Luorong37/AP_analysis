function tests = test_aaa_remaining_runtime
tests = functiontests(localfunctions);
end
function setupOnce(testCase)
unit_dir=fileparts(mfilename('fullpath')); app_dir=fileparts(fileparts(unit_dir));
old_path=path; restoredefaultpath; addpath(app_dir); AAA_startup();
testCase.TestData.oldPath=old_path;
end
function teardownOnce(testCase), path(testCase.TestData.oldPath); end

function testComparisonPreservesIndependentTimeAxes(testCase)
v=results_with([1:10;2:11]',[2:11;3:12]');
c=results_with([1:6;2:7]',[2:7;3:8]');
vp=profile("voltage",10); cp=profile("calcium",5);
r=compare_channel_traces(v,c,vp,cp,struct( ...
    'voltage_polarity',-1,'calcium_polarity',1));
verifyEqual(testCase,numel(r.sensitivity.voltage_time),10);
verifyEqual(testCase,numel(r.sensitivity.calcium_time),6);
verifyEqual(testCase,r.sensitivity.voltage_display,-r.sensitivity.voltage_metric);
end

function testFrequencyFindsSyntheticPeak(testCase)
fs=100; t=(0:999)'/fs; trace=sin(2*pi*7*t);
spec=struct('min_freq_hz',.5,'max_freq_hz',30, ...
    'wavelet_voices_per_octave',12,'wavelet_name',"amor");
r=analyze_frequency(trace,profile("voltage",fs),spec);
verifyEqual(testCase,r.roi.peak_frequency_hz,7,'AbsTol',.11);
end

function testStimMetricsMatchFrozenMeans(testCase)
trace=(1:20)'; windows=struct('baseline_frames',[1 5;11 15], ...
    'stim_frames',[6 10;16 20]);
r=analyze_stimulus_trials(trace,profile("voltage",10),windows,1);
verifyEqual(testCase,r.delta_mean,[5 5]);
verifyEqual(testCase,r.delta_auc,[2 2]);
end

function testRunSingleAllSkipWritesManifest(testCase)
output=string(tempname); cleanup=onCleanup(@() remove_tree(output)); %#ok<NASGU>
request=aaa.schema.default_request("dual","single");
request.input_path=string(tempdir); request.workflow.output_path=output;
request.workflow.preset="custom";
catalog=aaa.schema.section_catalog("dual");
sections=struct();
for name=[catalog.name], sections.(char(name))="skip"; end
request.workflow.sections=sections;
receipt=aaa.workflows.run_single(request,struct());
verifyEqual(testCase,receipt.status,"completed");
verifyTrue(testCase,isfile(fullfile(output,'analysis_manifest.mat')));
verifyEqual(testCase,numel(receipt.skipped_sections),numel(catalog));
end

function value=results_with(sensitivity,snr)
value=struct('trace_results',struct( ...
    'sensitivity',struct('data',sensitivity), ...
    'snr',struct('data',snr)));
end
function value=profile(role,fs)
value=struct('role',role,'frame_rate',fs);
end
function remove_tree(value)
if isfolder(value),rmdir(value,'s');end
end

function tests = test_aaa_roi_section
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
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testDualRunPersistsLegacyFilesAndRawParents(testCase)
output = string(tempname); cleanup = onCleanup(@() remove_tree(output)); %#ok<NASGU>
request = aaa.schema.default_request("dual","single");
request.input_path = string(tempdir);
request.workflow.output_path = output;
request.params.dual.map_bin = 1;
runtime = struct('roi_spec',struct('polygons',struct( ...
    'role',"voltage",'position',[3 3; 6 3; 6 6; 3 6])));
ctx = synthetic_context(request,runtime,[8 9],[5 7]);
[ctx,outcome] = aaa.sections.common.roi(ctx,"run");
[ctx,~]=save_analysis_results(ctx,"roi",outcome);
verifyEqual(testCase,outcome.state,"completed");
verifyTrue(testCase,isfile(fullfile(output,'results','roi.mat')));
for idx = 1:2
    raw = ctx.channels(idx).results.trace_results.raw;
    verifyEmpty(testCase,raw.info.parent_results);
    verifyEqual(testCase,size(raw.data,2),1);
    expected_frames = [5 7];
    verifyEqual(testCase,size(raw.data,1),expected_frames(idx));
end
end

function testApRunUsesSameSelectorAndOutputContract(testCase)
output = string(tempname); cleanup = onCleanup(@() remove_tree(output)); %#ok<NASGU>
request = aaa.schema.default_request("ap","single");
request.input_path = string(tempdir);
request.workflow.output_path = output;
request.params.ap.bin = 1;
runtime = struct('roi_spec',struct('polygons',struct( ...
    'role',"voltage",'position',[2 2; 5 2; 5 5; 2 5])));
ctx = synthetic_context(request,runtime,[7 8],6);
[ctx,outcome] = aaa.sections.common.roi(ctx,"run");
[ctx,~]=save_analysis_results(ctx,"roi",outcome);
verifyTrue(testCase,isfile(fullfile(output,'results','roi.mat')));
verifyEqual(testCase,size(ctx.channels(1).results.trace_results.raw.data),[6 1]);
end

function testDualReuseSkipsRedundantMapComputation(testCase)
source = string(tempname); output = string(tempname);
mkdir(source); cleanup = onCleanup(@() remove_both(source,output)); %#ok<NASGU>
rois = struct('bwmask',zeros(8,9),'bwmask_ca',zeros(8,9));
rois.bwmask(2:4,3:5) = 1; rois.bwmask_ca(2:4,3:5) = 1;
offset = [0 0]; %#ok<NASGU>
save(fullfile(source,'1_dual_roi_results.mat'),'rois','offset');
request = aaa.schema.default_request("dual","single");
request.input_path = string(tempdir);
request.workflow.output_path = output;
request.workflow.roi_path = source;
request.workflow.sections.registration = "skip";
request.workflow.sections.roi = "reuse";
ctx = synthetic_context(request,struct(),[8 9],[5 7]);
[ctx,outcome] = aaa.sections.common.roi(ctx,"reuse");
[ctx,~]=save_analysis_results(ctx,"roi",outcome);
verifyEqual(testCase,outcome.state,"completed");
verifyFalse(testCase,isfile(fullfile(output,'0_dual_sensitivity_map.mat')));
verifyTrue(testCase,isfile(fullfile(output,'results','roi.mat')));
end

function ctx = synthetic_context(request,runtime,frame_size,frame_counts)
ctx = aaa.sections.create_context(request,runtime);
for idx = 1:numel(ctx.channels)
    frames = frame_counts(min(idx,numel(frame_counts)));
    movie = reshape(single(1:(prod(frame_size)*frames)), ...
        frame_size(1),frame_size(2),frames);
    role = string(ctx.channels(idx).profile.role);
    ctx.channels(idx).data.movie_3d = movie;
    ctx.channels(idx).profile.roi.frame_size = frame_size;
    movie_info = struct('role',role,'frame_size',frame_size, ...
        'analysis_frame_size',frame_size,'frame_count',frames, ...
        'frame_rate',ctx.channels(idx).profile.frame_rate, ...
        'source_path',"synthetic", ...
        'transpose_before_analysis',false, ...
        'motion',struct('applied',false));
    ctx.channels(idx).data.movie_info = movie_info;
    ctx.channels(idx).results.movie_info = movie_info;
    ctx.channels(idx).results.trace_results = struct();
end
end

function remove_both(a,b)
remove_tree(a); remove_tree(b);
end

function remove_tree(path_value)
if isfolder(path_value), rmdir(path_value,'s'); end
end

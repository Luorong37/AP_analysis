function tests = test_aaa_select_rois
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

function testSingleMaskReplayUsesProfileRole(testCase)
movie = reshape(1:(6*7*9),6,7,9);
mask = zeros(6,7); mask(2:3,3:5) = 1;
[rois,traces,info] = select_rois({movie},profile("voltage",[6 7]), ...
    struct('masks',struct('voltage',mask)));
expected = squeeze(mean(movie(2:3,3:5,:),[1 2]))';
verifyEqual(testCase,traces.voltage,expected(:));
verifyEqual(testCase,rois.bwmask,mask);
verifyEqual(testCase,info.selection_mode,"reused_masks");
end

function testDualPolygonMapsByOffsetWithoutFrameAlignment(testCase)
voltage_movie = reshape(1:(8*9*5),8,9,5);
calcium_movie = reshape(1001:(1000+8*9*7),8,9,7);
profiles = [profile("voltage",[8 9]),profile("calcium",[8 9])];
polygon = struct('role',"voltage", ...
    'position',[4 3; 6 3; 6 5; 4 5]);
[rois,traces,info] = select_rois( ...
    {voltage_movie,calcium_movie},profiles, ...
    struct('polygons',polygon,'offset_xy',[2 1]));
verifyEqual(testCase,size(traces.voltage,1),5);
verifyEqual(testCase,size(traces.calcium,1),7);
verifyEqual(testCase,rois.bwmask_ca, ...
    imtranslate(rois.bwmask,[-2 -1],'nearest','FillValues',0));
verifyEqual(testCase,info.trace_alignment_rule, ...
    "none; each channel keeps its own frame count");
end

function testOneDualMaskBuildsPartnerMask(testCase)
movie = zeros(8,9,4);
voltage_mask = zeros(8,9); voltage_mask(3:5,4:6) = 1;
profiles = [profile("voltage",[8 9]),profile("calcium",[8 9])];
[rois,~,~] = select_rois({movie,movie},profiles, ...
    struct('masks',struct('voltage',voltage_mask),'offset_xy',[1 -2]));
expected = imtranslate(voltage_mask,[-1 2],'nearest','FillValues',0);
verifyEqual(testCase,rois.bwmask_ca,expected);
end

function value = profile(role,frame_size)
value = struct('role',role,'frame_rate',100, ...
    'transpose_before_analysis',false, ...
    'roi',struct('frame_size',frame_size));
end

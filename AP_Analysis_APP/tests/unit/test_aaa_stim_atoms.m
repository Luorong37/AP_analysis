function tests=test_aaa_stim_atoms
tests=functiontests(localfunctions);
end
function setupOnce(testCase)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(fullfile(root,'backend'),fullfile(root,'functions'));testCase.TestData.root=root;
end
function teardownOnce(testCase)
rmpath(fullfile(testCase.TestData.root,'backend'),fullfile(testCase.TestData.root,'functions'));
end
function testTrialMetricsSupportPerRoiPolarity(testCase)
trace=[1 4;2 3;4 2;6 1];profile=struct('role',"voltage",'frame_rate',2);w=struct('baseline_frames',[1 2],'stim_frames',[3 4]);
m=analyze_stimulus_trials(trace,profile,w,[1 -1]);
verifyEqual(testCase,m.delta_mean,[3.5;2],'AbsTol',1e-12);
end
function testPeakCounting(testCase)
w=struct('baseline_frames',[1 3;5 6],'stim_frames',[4 4;7 10]);r=count_stimulus_peaks({[2 4 8],[1 9]},w);
verifyEqual(testCase,r.baseline,[1 0;1 0]);verifyEqual(testCase,r.stim,[1 1;0 1]);
end
function testGratingTuning(testCase)
stim=[5 1 4 2;2 3 1 4];base=zeros(2,4);angles=[0 180 90 270];t=analyze_grating_tuning(stim,base,angles);
verifyEqual(testCase,t.pref_dir,[0;270]);verifySize(testCase,t.response_by_condition,[2 4]);verifyTrue(testCase,all(isfinite(t.gdsi)));
end
function testAutomaticSingleChannelDiscovery(testCase)
folder=tempname;mkdir(folder);cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
method_file=fullfile(folder,'method_manifest.mat');logs_file=fullfile(folder,'logs.mat');
manifest=struct('spec',struct('stimSpec',struct('selectedProgram',"drifting grating",'selectedLabel',"grating",'isi',.2,'duration',.2,'orientations',[0 90])), ...
    'actual',struct('stimRuntime',struct('ifi',.1)));save(method_file,'manifest');
manifest=struct('spec',struct('recordmode',"visualstim"),'refs',struct('method_manifest',method_file),'artifacts',struct('logs_mat',logs_file));save(fullfile(folder,'cycle_manifest.mat'),'manifest');
n=12;logs=struct();logs.sync=table((0:n-1)'*.1,zeros(n,1),zeros(n,1),(1:n)',(1:n)', ...
    'VariableNames',{'PTB_VBL_Time','Unused1','Unused2','Camera1_Frame','Camera2_Frame'});save(logs_file,'logs');
channel=struct('camera_index',1,'profile',struct('role',"voltage",'frame_rate',10), ...
    'data',struct('movie_info',struct('frame_count',20)),'results',struct());
[context,w]=resolve_visualstim(folder,channel,struct(),struct());
verifyTrue(testCase,context.supported);verifyTrue(testCase,w.is_grating);verifyEqual(testCase,w.voltage.baseline_frames,[1 2;5 6]);verifyEqual(testCase,w.voltage.stim_frames,[3 4;7 8]);
end

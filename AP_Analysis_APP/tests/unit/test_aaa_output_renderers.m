function tests=test_aaa_output_renderers
tests=functiontests(localfunctions);
end
function setupOnce(testCase)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(fullfile(root,'backend'),fullfile(root,'functions'));testCase.TestData.root=root;
end
function teardownOnce(testCase),rmpath(fullfile(testCase.TestData.root,'backend'),fullfile(testCase.TestData.root,'functions'));end
function testDualOverviewOutputContract(testCase)
[folder,cleanup]=new_folder(); %#ok<ASGLU>
[channels,maps,rois]=fixture();
bv=struct('background_raw',rand(20,1));bv.background_fitted=bv.background_raw;
bc=struct('background_raw',rand(20,1));bc.background_fitted=bc.background_raw;
shared=struct('maps',maps,'roi',struct('rois',rois), ...
    'trace_details',struct('background',struct('voltage',bv,'calcium',bc)));
render_analysis_overview(channels,folder,struct('mode',"dual",'voltage_polarity',-1,'calcium_polarity',1),shared);
expected=["0_dual_sensitivity_map.png","0_dual_sensitivity_map_with_roi.png","0_dual_average_color_merge.tif","0_dual_average_color_merge_with_roi.tif","1_dual_raw_trace.png","1_dual_background_correction_summary.png","2_dual_bleach_correction_stacked.png","3_dual_noise_reference_comparison.png","4_dual_trace_summary.png"];
verifyTrue(testCase,all(isfile(fullfile(folder,expected))));clear cleanup
end
function testComparisonAndStimOutputContract(testCase)
[folder,cleanup]=new_folder();[channels,~,~]=fixture(); %#ok<ASGLU>
comparison=compare_channel_traces(channels(1).results,channels(2).results,channels(1).profile,channels(2).profile,struct('voltage_polarity',-1,'calcium_polarity',1));render_dual_comparison(comparison,folder,struct());
verifyTrue(testCase,isfile(fullfile(folder,'5_dual_sensitivity_overlap.png')));verifyTrue(testCase,isfile(fullfile(folder,'5_dual_snr_overlap_rois','ROI_001.fig')));
w=struct('supported',true,'is_grating',true,'orientations',[0;90],'condition_index',[1;2],'condition_labels',["0";"90"],'voltage',struct('stim_time_ranges',[.5 .8;1.2 1.5]));w.calcium=w.voltage;
t=analyze_grating_tuning([1 2],[0 0],[0 90]);stim=struct('supported',true,'windows',w,'tuning',struct('voltage',t,'calcium',t));render_visualstim_results(stim,channels,folder,struct('mode',"dual"));
verifyTrue(testCase,isfile(fullfile(folder,'4_dual_trace_summary_with_stim.png')));verifyTrue(testCase,isfile(fullfile(folder,'7_stim_tuning_by_roi','roi_001_stim_tuning.fig')));clear cleanup
end
function testFrequencyOutputContract(testCase)
[folder,cleanup]=new_folder();fs=20;t=(0:127)'/fs;spec=struct('min_freq_hz',.5,'max_freq_hz',9,'wavelet_voices_per_octave',12,'wavelet_name',"amor");
v=analyze_frequency(sin(2*pi*3*t),struct('role',"voltage",'frame_rate',fs),spec);c=analyze_frequency(sin(2*pi*2*t),struct('role',"calcium",'frame_rate',fs),spec);results=struct('voltage',v,'calcium',c);render_frequency_results(results,folder,struct());
verifyTrue(testCase,isfile(fullfile(folder,'8_fourier_summary.png')));verifyTrue(testCase,isfile(fullfile(folder,'8_wavelet_summary.fig')));verifyTrue(testCase,isfile(fullfile(folder,'8_wavelet_roi_pairs','ROI_001.png')));clear cleanup
end
function [channels,maps,rois]=fixture()
rng(1);movie=rand(8,8,20);mask=zeros(8);mask(3:5,3:5)=1;roles=["voltage","calcium"];channels=repmat(struct('profile',struct(),'data',struct(),'results',struct()),1,2);maps=struct();
for i=1:2,role=roles(i);x=randn(20,1);tr=struct();for name=["raw","bleach_removed","noise_reference","sensitivity","snr","raw_smoothed","sensitivity_smoothed","snr_smoothed"],tr.(name)=struct('data',x);end;channels(i)=struct('profile',struct('role',role,'frame_rate',10),'data',struct('movie_3d',movie,'roi_mask',mask),'results',struct('trace_results',tr));maps.(role)=randn(8);end
rois=struct('bwmask',mask,'bwmask_ca',mask);
end
function [folder,cleanup]=new_folder(),folder=tempname;mkdir(folder);cleanup=onCleanup(@()rmdir(folder,'s'));end

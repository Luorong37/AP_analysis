function response = count_stimulus_peaks(peak_indices, windows)
%COUNT_STIMULUS_PEAKS Count accepted ROI peaks in baseline/stim ranges.
% Inputs:
%   peak_indices 1-by-ROI cell array of one-based accepted frame indices.
%   windows.baseline_frames/windows.stim_frames Matching [trial x 2]
%       inclusive one-based frame ranges.
% Output:
%   response Baseline, stimulus, and delta counts as [ROI x trial].
if ~iscell(peak_indices) || ~isstruct(windows) ...
        || ~all(isfield(windows,{'baseline_frames','stim_frames'}))
    error('AAA:Functions:InvalidPeakStimInput', ...
        'Peak indices and baseline/stim windows are required.');
end
nrois=numel(peak_indices);ntrials=size(windows.stim_frames,1);
baseline=zeros(nrois,ntrials);stim=zeros(nrois,ntrials);
for roi=1:nrois
    indices=double(peak_indices{roi}(:));
    for trial=1:ntrials
        b=windows.baseline_frames(trial,:);s=windows.stim_frames(trial,:);
        baseline(roi,trial)=sum(indices>=b(1)&indices<=b(2));
        stim(roi,trial)=sum(indices>=s(1)&indices<=s(2));
    end
end
response=struct('baseline',baseline,'stim',stim,'delta',stim-baseline, ...
    'units',"accepted peak count",'windows',windows);
end

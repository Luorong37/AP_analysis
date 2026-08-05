function metrics = analyze_stimulus_trials(trace_matrix,profile,windows,polarity)
%ANALYZE_STIMULUS_TRIALS Frozen Dual baseline-vs-stim trial metrics.
% Inputs:
%   trace_matrix Numeric [frames x ROI] matrix.
%   profile Scalar bound profile requiring role and frame_rate in Hz.
%   windows.baseline_frames and windows.stim_frames are matching
%       [trial x 2] inclusive one-based integer frame ranges.
%   polarity Nonzero scalar or one value per ROI; values normalize to +/-1.
% Output:
%   metrics Per-ROI/per-trial values stored as [ROI x trial], plus summary,
%           effective polarity, windows, role, and frame rate.
required={'baseline_frames','stim_frames'};
if ~isstruct(windows) || ~all(isfield(windows,required))
    error('AAA:Functions:InvalidStimWindows', ...
        'Stim windows require baseline_frames and stim_frames.');
end
baseline=double(windows.baseline_frames); stimulus=double(windows.stim_frames);
if size(baseline,2)~=2 || ~isequal(size(baseline),size(stimulus)) ...
        || any(~isfinite([baseline(:);stimulus(:)])) ...
        || any([baseline(:);stimulus(:)]~=round([baseline(:);stimulus(:)])) ...
        || any([baseline(:);stimulus(:)]<1) ...
        || any([baseline(:);stimulus(:)]>size(trace_matrix,1)) ...
        || any(baseline(:,1)>baseline(:,2)) || any(stimulus(:,1)>stimulus(:,2))
    error('AAA:Functions:InvalidStimWindows', ...
        'Stim frame ranges must be valid inclusive integer ranges.');
end
trace=double(trace_matrix);nrois=size(trace,2);polarity=normalize_polarity(polarity,nrois);trace=trace.*reshape(polarity,1,[]);ntrials=size(stimulus,1);
names={'baseline_mean','stim_mean','delta_mean','baseline_peak','stim_peak', ...
    'delta_peak','baseline_auc','stim_auc','delta_auc'};
metrics=struct('role',string(profile.role),'frame_rate',double(profile.frame_rate), ...
    'polarity',double(polarity));
for idx=1:numel(names), metrics.(names{idx})=NaN(nrois,ntrials); end
for trial=1:ntrials
    b=trace(baseline(trial,1):baseline(trial,2),:);
    s=trace(stimulus(trial,1):stimulus(trial,2),:);
    metrics.baseline_mean(:,trial)=mean(b,1,'omitnan')';
    metrics.stim_mean(:,trial)=mean(s,1,'omitnan')';
    metrics.delta_mean(:,trial)=metrics.stim_mean(:,trial)-metrics.baseline_mean(:,trial);
    metrics.baseline_peak(:,trial)=max(b,[],1)';
    metrics.stim_peak(:,trial)=max(s,[],1)';
    metrics.delta_peak(:,trial)=metrics.stim_peak(:,trial)-metrics.baseline_peak(:,trial);
    metrics.baseline_auc(:,trial)=trapz(b,1)'/double(profile.frame_rate);
    metrics.stim_auc(:,trial)=trapz(s,1)'/double(profile.frame_rate);
    metrics.delta_auc(:,trial)=metrics.stim_auc(:,trial)-metrics.baseline_auc(:,trial);
end
metrics.summary=repmat(struct(),nrois,1);
for roi=1:nrois
    metrics.summary(roi).delta_mean=mean(metrics.delta_mean(roi,:),'omitnan');
    metrics.summary(roi).delta_peak=mean(metrics.delta_peak(roi,:),'omitnan');
    metrics.summary(roi).delta_auc=mean(metrics.delta_auc(roi,:),'omitnan');
    metrics.summary(roi).p_mean=paired_p(metrics.baseline_mean(roi,:),metrics.stim_mean(roi,:));
    metrics.summary(roi).p_peak=paired_p(metrics.baseline_peak(roi,:),metrics.stim_peak(roi,:));
    metrics.summary(roi).p_auc=paired_p(metrics.baseline_auc(roi,:),metrics.stim_auc(roi,:));
end
metrics.windows=windows;
end
function value=normalize_polarity(value,nrois)
value=double(value(:)');if isscalar(value),value=repmat(value,1,nrois);end
if numel(value)~=nrois||any(~isfinite(value))||any(value==0),error('AAA:Functions:InvalidStimPolarity','Polarity must be scalar or one nonzero value per ROI.');end
value(value>0)=1;value(value<0)=-1;
end
function p=paired_p(x,y)
valid=isfinite(x)&isfinite(y); x=x(valid); y=y(valid);
if numel(x)<2, p=NaN; return; end
try, p=signrank(x,y); catch, [~,p]=ttest(x,y); end
end

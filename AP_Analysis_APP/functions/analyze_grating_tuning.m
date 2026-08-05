function tuning = analyze_grating_tuning(stim_response, baseline_response, orientations)
%ANALYZE_GRATING_TUNING Compute frozen Dual direction/orientation metrics.
% Inputs:
%   stim_response/baseline_response Numeric [ROI x trial] matrices with
%       identical size.
%   orientations Numeric [trial x 1] directions in degrees.
% Output:
%   tuning Per-condition response matrices and per-ROI preferred direction,
%          preferred orientation, DSI/OSI, and global DSI/OSI.
stim_response=double(stim_response);baseline_response=double(baseline_response);
orientations=double(orientations(:));
if size(stim_response,2)~=numel(orientations) ...
        || ~isequal(size(stim_response),size(baseline_response))
    error('AAA:Functions:InvalidGratingResponse', ...
        'Responses must be ROI-by-trial matrices matching orientations.');
end
unique_angles=unique(orientations,'sorted');nrois=size(stim_response,1);
tuning=struct('angles',unique_angles(:)', ...
    'trial_response',max(stim_response,0), ...
    'nonstim_trial_response',max(baseline_response,0), ...
    'response_by_condition',NaN(nrois,numel(unique_angles)), ...
    'nonstim_response_by_condition',NaN(nrois,numel(unique_angles)), ...
    'pref_dir',NaN(nrois,1),'pref_ori',NaN(nrois,1), ...
    'dsi',NaN(nrois,1),'osi',NaN(nrois,1), ...
    'gdsi',NaN(nrois,1),'gosi',NaN(nrois,1));
for cond=1:numel(unique_angles)
    mask=orientations==unique_angles(cond);
    tuning.response_by_condition(:,cond)=mean(tuning.trial_response(:,mask),2,'omitnan');
    tuning.nonstim_response_by_condition(:,cond)=mean(tuning.nonstim_trial_response(:,mask),2,'omitnan');
end
theta=deg2rad(unique_angles(:)');
for roi=1:nrois
    r=tuning.response_by_condition(roi,:);valid=isfinite(r);
    if ~any(valid),continue;end
    [~,pref]=max(r);tuning.pref_dir(roi)=unique_angles(pref);tuning.pref_ori(roi)=mod(unique_angles(pref),180);
    denom=max(eps,sum(r(valid)));
    tuning.gdsi(roi)=abs(sum(r(valid).*exp(1i*theta(valid))))/denom;
    tuning.gosi(roi)=abs(sum(r(valid).*exp(2i*theta(valid))))/denom;
    opp=closest_angle(unique_angles,mod(unique_angles(pref)+180,360),360);
    ortho=closest_angle(unique_angles,mod(unique_angles(pref)+90,180),180);
    tuning.dsi(roi)=(r(pref)-r(opp))/max(eps,r(pref)+r(opp));
    tuning.osi(roi)=(r(pref)-r(ortho))/max(eps,r(pref)+r(ortho));
end
end
function idx=closest_angle(values,target,period)
d=abs(mod(double(values)-target+period/2,period)-period/2);[~,idx]=min(d);
end

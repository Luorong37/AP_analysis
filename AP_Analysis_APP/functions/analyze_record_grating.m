function summary = analyze_record_grating(entries,roles)
%ANALYZE_RECORD_GRATING Aggregate saved Cycle tuning without mixing units.
% Curves are aligned by the union of directions modulo 360. Missing
% directions remain NaN. Per-Cycle metrics and metrics recomputed from the
% Cycle-mean curve are both retained.

roles = lower(string(roles(:)));
summary = struct('status',"not_available", ...
    'reason',"No compatible grating tuning results were found.", ...
    'included_cycles',strings(0,1),'excluded_cycles',strings(0,1), ...
    'orientations',[],'channels',struct(), ...
    'voltage',struct(),'voltage_peak_count',struct(), ...
    'voltage_sensitivity',struct(),'calcium',struct(), ...
    'response_method_by_cycle',struct(), ...
    'curve_rule',struct( ...
        'orientation_rule',"union after modulo 360; missing directions remain NaN", ...
        'cycle_average_rule',"mean across available cycles with omitnan", ...
        'metric_mean_rule',"mean and SEM of per-cycle metrics", ...
        'metric_from_mean_curve_rule',"metrics recomputed from the cycle-mean grating curve", ...
        'nonstim_role',"saved and plotted but excluded from DSI/OSI", ...
        'unit_rule',"peak count and mean sensitivity are never pooled"));
if isempty(entries), return; end

all_angles = [];
valid_any = false(numel(entries),1);
for idx = 1:numel(entries)
    tuning = entries(idx).grating_tuning;
    names = fieldnames(tuning);
    for j = 1:numel(names)
        candidate = tuning.(names{j});
        if complete_tuning(candidate)
            all_angles = [all_angles,tuning_angles(candidate)]; %#ok<AGROW>
            valid_any(idx) = true;
        end
    end
end
if isempty(all_angles), return; end
angles = sort(unique(mod(double(all_angles(:)'),360)));
summary.orientations = angles;
summary.included_cycles = string({entries(valid_any).cycle_name})';
summary.excluded_cycles = string({entries(~valid_any).cycle_name})';

for role = roles'
    role_name = char(role);
    if role == "voltage"
        [peak_mask,sensitivity_mask,methods] = voltage_method_masks(entries);
        summary.response_method_by_cycle.voltage = table( ...
            string({entries.cycle_name})',methods, ...
            'VariableNames',{'Cycle','EffectiveMethod'});
        main_mask = reshape(arrayfun(@(entry)has_tuning(entry,'voltage'),entries),[],1);
        if any(main_mask)
            main_methods = unique(methods(main_mask));
            if numel(main_methods)==1
                summary.voltage = average_variant(entries,main_mask,'voltage',angles);
            else
                summary.voltage = struct('status',"mixed_units", ...
                    'reason',"Primary voltage tuning contains multiple effective response methods; use voltage_peak_count or voltage_sensitivity.", ...
                    'methods',main_methods);
            end
        end
        if any(peak_mask)
            summary.voltage_peak_count = average_variant(entries,peak_mask,'voltage',angles);
        else
            summary.voltage_peak_count = unavailable("No peak-count tuning cycles were available.");
        end
        explicit_sensitivity = reshape(arrayfun( ...
            @(entry)has_tuning(entry,'voltage_sensitivity'),entries),[],1);
        sensitivity_mask = sensitivity_mask | explicit_sensitivity;
        if any(sensitivity_mask)
            field_by_cycle = repmat("voltage",numel(entries),1);
            field_by_cycle(explicit_sensitivity) = "voltage_sensitivity";
            summary.voltage_sensitivity = average_variant_per_cycle( ...
                entries,sensitivity_mask,field_by_cycle,angles);
        else
            summary.voltage_sensitivity = unavailable("No mean-sensitivity tuning cycles were available.");
        end
        summary.channels.voltage = struct( ...
            'primary',summary.voltage,'peak_count',summary.voltage_peak_count, ...
            'sensitivity',summary.voltage_sensitivity);
    else
        mask = reshape(arrayfun(@(entry)has_tuning(entry,role_name),entries),[],1);
        if any(mask)
            value = average_variant(entries,mask,role_name,angles);
        else
            value = unavailable("No compatible "+role+" tuning cycles were available.");
        end
        summary.channels.(role_name) = value;
        summary.(role_name) = value;
    end
end
summary.status = "completed";
summary.reason = "";
end

function output = average_variant(entries,mask,field_name,angles)
fields = repmat(string(field_name),numel(entries),1);
output = average_variant_per_cycle(entries,mask,fields,angles);
end

function output = average_variant_per_cycle(entries,mask,fields,angles)
selected = find(mask(:));
first = entries(selected(1)).grating_tuning.(char(fields(selected(1))));
nrois = size(first.response_by_condition,1);
nangles = numel(angles); ncycles = numel(selected);
stim = NaN(nrois,nangles,ncycles);
nonstim = NaN(nrois,nangles,ncycles);
metrics = initialize_metrics(nrois,ncycles);
included = strings(ncycles,1);
for out_idx = 1:ncycles
    idx = selected(out_idx); included(out_idx)=entries(idx).cycle_name;
    tuning = entries(idx).grating_tuning.(char(fields(idx)));
    if size(tuning.response_by_condition,1) ~= nrois
        error('AAA:Record:TuningROIMismatch', ...
            'ROI count mismatch in grating tuning for cycle %s.',entries(idx).cycle_name);
    end
    source_angles = tuning_angles(tuning);
    for angle_idx = 1:numel(source_angles)
        target = find(abs(angles-source_angles(angle_idx))<1e-9,1);
        if isempty(target), continue; end
        stim(:,target,out_idx) = tuning.response_by_condition(:,angle_idx);
        nonstim(:,target,out_idx) = tuning.nonstim_response_by_condition(:,angle_idx);
    end
    names = fieldnames(metrics);
    for metric_idx = 1:numel(names)
        name = names{metric_idx};
        if isfield(tuning,name), metrics.(name)(:,out_idx)=tuning.(name)(:); end
    end
end
stim_mean = mean(stim,3,'omitnan'); nonstim_mean = mean(nonstim,3,'omitnan');
stim_n = sum(isfinite(stim),3); nonstim_n = sum(isfinite(nonstim),3);
output = struct('status',"completed",'included_cycles',included, ...
    'stim_per_cycle',stim,'nonstim_per_cycle',nonstim, ...
    'stim_mean',stim_mean,'stim_sem',std(stim,0,3,'omitnan')./sqrt(max(1,stim_n)), ...
    'stim_n',stim_n,'nonstim_mean',nonstim_mean, ...
    'nonstim_sem',std(nonstim,0,3,'omitnan')./sqrt(max(1,nonstim_n)), ...
    'nonstim_n',nonstim_n,'metrics_per_cycle',metrics, ...
    'metrics_mean_across_cycles',summarize_metrics(metrics), ...
    'metrics_from_mean_curve',compute_metrics(stim_mean,angles));
end

function [peak_mask,sensitivity_mask,methods] = voltage_method_masks(entries)
methods = strings(numel(entries),1);
for idx = 1:numel(entries)
    method = "";
    if isstruct(entries(idx).voltage_response) ...
            && isfield(entries(idx).voltage_response,'effective')
        method = lower(string(entries(idx).voltage_response.effective));
    elseif isstruct(entries(idx).voltage_response) ...
            && isfield(entries(idx).voltage_response,'units') ...
            && contains(lower(string(entries(idx).voltage_response.units)),"sensitivity")
        method = "mean_sensitivity";
    elseif has_tuning(entries(idx),'voltage')
        method = "peak_count_legacy";
    end
    methods(idx) = method;
end
peak_mask = reshape(arrayfun(@(entry)has_tuning(entry,'voltage'),entries),[],1) ...
    & ismember(methods,["peak_count","peak_count_legacy"]);
sensitivity_mask = reshape(arrayfun(@(entry)has_tuning(entry,'voltage'),entries),[],1) ...
    & methods=="mean_sensitivity";
end

function tf = has_tuning(entry,name)
tf = isstruct(entry.grating_tuning) && isfield(entry.grating_tuning,name) ...
    && complete_tuning(entry.grating_tuning.(name));
end
function tf = complete_tuning(value)
tf = isstruct(value) && isfield(value,'response_by_condition') ...
    && isfield(value,'nonstim_response_by_condition') ...
    && (isfield(value,'angles') || isfield(value,'unique_orientations'));
end
function angles = tuning_angles(value)
if isfield(value,'angles'), angles=value.angles;else,angles=value.unique_orientations;end
angles = mod(double(angles(:)'),360);
end
function metrics = initialize_metrics(nrois,ncycles)
names={'pref_dir','pref_ori','gdsi','gosi','dsi','osi'}; metrics=struct();
for idx=1:numel(names),metrics.(names{idx})=NaN(nrois,ncycles);end
end
function summary = summarize_metrics(metrics)
summary=struct(); names=fieldnames(metrics);
for idx=1:numel(names)
    values=metrics.(names{idx});
    summary.(names{idx})=struct('mean',mean(values,2,'omitnan'), ...
        'sem',std(values,0,2,'omitnan')./sqrt(max(1,sum(isfinite(values),2))), ...
        'n',sum(isfinite(values),2));
end
end
function metrics = compute_metrics(response,angles)
nrois=size(response,1); names={'pref_dir','pref_ori','gdsi','gosi','dsi','osi'};
metrics=struct();for idx=1:numel(names),metrics.(names{idx})=NaN(nrois,1);end
theta=deg2rad(angles(:)');
for roi=1:nrois
    r=response(roi,:);
    if all(~isfinite(r))||sum(r,'omitnan')<=0,continue;end
    [~,pref]=max(r,[],'omitnan');metrics.pref_dir(roi)=angles(pref);metrics.pref_ori(roi)=mod(angles(pref),180);
    denom=sum(r,'omitnan');metrics.gdsi(roi)=abs(sum(r.*exp(1i*theta),'omitnan')/denom);metrics.gosi(roi)=abs(sum(r.*exp(2i*theta),'omitnan')/denom);
    opp=closest(angles,mod(angles(pref)+180,360));ortho=closest(angles,mod(angles(pref)+90,360));
    metrics.dsi(roi)=(r(pref)-r(opp))/max(eps,r(pref)+r(opp));
    metrics.osi(roi)=(r(pref)-r(ortho))/max(eps,r(pref)+r(ortho));
end
end
function idx=closest(angles,target)
[~,idx]=min(abs(mod(angles-target+180,360)-180));
end
function value=unavailable(reason)
value=struct('status',"not_available",'reason',string(reason));
end

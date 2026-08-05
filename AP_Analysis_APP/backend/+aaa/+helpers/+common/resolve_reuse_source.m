function [source_path, report] = resolve_reuse_source(mode, search_root, plan, target, options)
%RESOLVE_REUSE_SOURCE Find the newest semantically compatible result.
% Filenames and manifests only locate candidate folders. Compatibility is
% decided from the saved result schema. Search is bounded to keep App path
% callbacks responsive and never reads movie arrays.

if nargin < 4 || isempty(target), target = struct(); end
if nargin < 5 || isempty(options), options = struct(); end
options = normalize_options(options);
mode = lower(strtrim(string(mode)));
search_root = strtrim(string(search_root));
if ~ismember(mode,["ap","dual"])
    error('AAA:Reuse:InvalidMode','mode must be ap or dual.');
end
if ~isscalar(search_root) || strlength(search_root)==0 || ~isfolder(search_root)
    source_path = "";
    report = empty_report(search_root,target,options);
    report.status = "search_root_missing";
    return;
end

report = empty_report(search_root,target,options);
[direct_ok,direct_reason,direct_inspection] = ...
    candidate_complete(search_root,mode,plan,options.discovery_max_mat_files);
if direct_ok && target_matches(search_root,target)
    source_path = search_root;
    report.status = "resolved_direct";
    report.chosen_path = source_path;
    report.selection_rule = "explicit folder passed semantic validation";
    report.chosen_inspection = direct_inspection;
    report.candidates = candidate_record(search_root,candidate_time(search_root), ...
        true,"complete",direct_inspection);
    return;
end

candidate_paths = strings(0,1);
roots = targeted_scan_roots(search_root,target);
for idx = 1:numel(roots)
    [found,scan] = scan_candidate_folders(roots(idx),options);
    candidate_paths = [candidate_paths;found]; %#ok<AGROW>
    report.scanned_directory_count = report.scanned_directory_count + scan.directories;
    report.scanned_anchor_count = report.scanned_anchor_count + scan.anchors;
    if scan.truncated
        report.truncated = true;
        report.warnings(end+1,1) = "Reuse scan limit reached under " + roots(idx) + ...
            "; automatic selection is disabled for this incomplete scan."; %#ok<AGROW>
    end
end
candidate_paths = unique(candidate_paths,'stable');

records = repmat(candidate_record("",NaT,false,"",struct()),0,1);
for idx = 1:numel(candidate_paths)
    candidate = candidate_paths(idx);
    if ~target_matches(candidate,target)
        records(end+1,1) = candidate_record(candidate,candidate_time(candidate), ...
            false,"different Cycle",struct()); %#ok<AGROW>
        continue;
    end
    [ok,reason,inspection] = candidate_complete(candidate,mode,plan, ...
        options.discovery_max_mat_files);
    records(end+1,1) = candidate_record(candidate,candidate_time(candidate), ...
        ok,reason,inspection); %#ok<AGROW>
end
report.candidates = records;
report.warnings = unique(report.warnings,'stable');
if report.truncated
    source_path = "";
    report.status = "scan_truncated";
    report.selection_rule = "incomplete bounded scan cannot auto-select";
    return;
end

valid_idx = find([records.valid]);
if isempty(valid_idx)
    source_path = "";
    report.status = "not_found";
    if strlength(direct_reason)>0
        report.selection_rule = "direct folder rejected: " + direct_reason;
    end
    return;
end
times = [records(valid_idx).timestamp];
newest_time=max(times);
ties=valid_idx(times==newest_time);
tie_paths=lower(string({records(ties).path}));
[~,order]=sort(tie_paths);
chosen=ties(order(1));
source_path = records(chosen).path;
report.status = "resolved_latest";
report.chosen_path = source_path;
report.chosen_inspection = records(chosen).inspection;
report.selection_rule = ...
    "newest semantically complete candidate by manifest/file time";
end

function options = normalize_options(options)
defaults = struct('discovery_max_depth',4,'discovery_max_mat_files',500);
names = fieldnames(defaults);
for idx = 1:numel(names)
    name = names{idx};
    if ~isfield(options,name) || isempty(options.(name))
        options.(name) = defaults.(name);
    end
    value = double(options.(name));
    if ~isscalar(value) || ~isfinite(value) || value < 1 || value ~= round(value)
        error('AAA:Reuse:InvalidDiscoveryLimit','%s must be a positive integer.',name);
    end
    options.(name) = value;
end
end

function [folders,report] = scan_candidate_folders(root,options)
folders = strings(0,1);
queue = struct('path',string(root),'depth',0);
head = 1; directories = 0; anchors = 0; truncated = false;
while head <= numel(queue)
    current = queue(head); head = head + 1;
    directories = directories + 1;
    listing = dir(current.path);
    files = listing(~[listing.isdir]);
    file_names = string({files.name});
    is_anchor = file_names == "analysis_manifest.mat" ...
        | endsWith(file_names,"_results.mat",'IgnoreCase',true);
    anchors = anchors + nnz(is_anchor);
    if any(is_anchor), folders(end+1,1) = current.path; end %#ok<AGROW>
    if anchors >= options.discovery_max_mat_files
        truncated = true; break;
    end
    if current.depth >= options.discovery_max_depth, continue; end
    children = listing([listing.isdir]);
    child_names = string({children.name});
    children = children(child_names~="." & child_names~=".." ...
        & ~startsWith(child_names,"."));
    for idx = 1:numel(children)
        queue(end+1) = struct('path', ...
            string(fullfile(children(idx).folder,children(idx).name)), ...
            'depth',current.depth+1); %#ok<AGROW>
    end
end
report = struct('directories',directories,'anchors',anchors,'truncated',truncated);
end

function roots = targeted_scan_roots(search_root,target)
roots = strings(0,1);
if target_matches(search_root,target), roots = search_root; return; end
cycle_name = ""; label = "";
if isfield(target,'cycle_name'), cycle_name=string(target.cycle_name); end
if isfield(target,'label'), label=string(target.label); end
if strlength(label)>0
    candidate = fullfile(search_root,char(replace(label,"/",filesep)));
    if isfolder(candidate), roots(end+1,1)=string(candidate); end %#ok<AGROW>
end
if strlength(cycle_name)>0
    candidate = fullfile(search_root,cycle_name);
    if isfolder(candidate), roots(end+1,1)=string(candidate); end %#ok<AGROW>
    records = dir(fullfile(search_root,'Rec*'));
    records = records([records.isdir]);
    for idx=1:numel(records)
        candidate = fullfile(records(idx).folder,records(idx).name,cycle_name);
        if isfolder(candidate), roots(end+1,1)=string(candidate); end %#ok<AGROW>
    end
end
roots = unique(roots,'stable');
if isempty(roots), roots=search_root; end
end

function [tf,reason,inspection] = candidate_complete(folder,mode,plan,max_files)
tf=false; reason=""; inspection=struct();
if ~isfolder(folder), reason="not a folder"; return; end
if ~aaa.io.manifest_allows_reuse(folder)
    reason="manifest status is not completed"; return;
end
if reused(plan,'input') || reused(plan,'trace')
    inspection = aaa.helpers.common.inspect_result_folder(folder,mode,plan,max_files);
    if ~inspection.valid, reason=inspection.reason; return; end
end
missing=strings(0,1);
semantic_names=["roi","registration","peak","stim", ...
    "comparison","frequency","visualization"];
needs_semantic=false;
for name=semantic_names,needs_semantic=needs_semantic||reused(plan,char(name));end
if needs_semantic
    try
        [results,~]=load_analysis_results(folder,struct( ...
            'mode',mode,'scope',"single",'max_mat_files',max_files));
    catch ME
        reason="result compatibility load failed: "+string(ME.message);return;
    end
else
    results=struct();
end
if reused(plan,'roi')&&~isfield(results,'roi'),missing(end+1,1)="roi";end %#ok<AGROW>
if reused(plan,'registration')&&~isfield(results,'registration'),missing(end+1,1)="registration";end %#ok<AGROW>
if reused(plan,'peak')&&(~isfield(results,'voltage') ...
        ||~isfield(results.voltage,'peak_results')),missing(end+1,1)="peak";end %#ok<AGROW>
if reused(plan,'stim')&&~isfield(results,'stim'),missing(end+1,1)="stim";end %#ok<AGROW>
if reused(plan,'comparison')&&~isfield(results,'comparison'),missing(end+1,1)="comparison";end %#ok<AGROW>
if reused(plan,'frequency')&&~isfield(results,'frequency'),missing(end+1,1)="frequency";end %#ok<AGROW>
if reused(plan,'visualization')&&~isfield(results,'visualization')
    missing(end+1,1)="visualization"; %#ok<AGROW>
end
if ~isempty(missing),reason="missing semantic result(s): "+strjoin(missing,", ");return;end
if reused(plan,'motion') && ~(isfile(fullfile(folder,'shared_motion_shifts_result.mat')) ...
        || isfile(fullfile(folder,'motion_shifts_result.mat')))
    reason="missing motion shift result"; return;
end
tf=true; reason="complete";
end

function tf = reused(plan,name)
tf=isstruct(plan)&&isfield(plan,name)&&string(plan.(name))=="reuse";
end

function tf = target_matches(candidate,target)
tf=true; names=strings(0,1);
if isfield(target,'cycle_name'), names(end+1,1)=string(target.cycle_name); end %#ok<AGROW>
if isfield(target,'label'), names(end+1,1)=string(target.label); end %#ok<AGROW>
names=unique(names(strlength(names)>0)); if isempty(names), return; end
parts=split(replace(string(candidate),'/','\'),'\');
tf=any(ismember(lower(parts),lower(names)));
end

function value = candidate_time(folder)
value=NaT; manifest_file=fullfile(folder,'analysis_manifest.mat');
if isfile(manifest_file)
    try
        manifest=aaa.io.load_manifest(manifest_file);
        if isfield(manifest,'completed_at') && ~isempty(manifest.completed_at)
            value=datetime(manifest.completed_at);
        end
    catch
    end
end
if isnat(value)
    listing=dir(folder);
    if ~isempty(listing), value=datetime(max([listing.datenum]),'ConvertFrom','datenum');
    else, value=datetime(0,'ConvertFrom','datenum'); end
end
end

function report = empty_report(search_root,target,options)
report=struct('status',"pending",'search_root',string(search_root), ...
    'target',target,'chosen_path',"",'selection_rule',"", ...
    'options',options,'truncated',false,'warnings',strings(0,1), ...
    'scanned_directory_count',0,'scanned_anchor_count',0, ...
    'chosen_inspection',struct(), ...
    'candidates',repmat(candidate_record("",NaT,false,"",struct()),0,1));
end

function record = candidate_record(path_value,timestamp,valid,reason,inspection)
record=struct('path',string(path_value),'timestamp',timestamp, ...
    'valid',logical(valid),'reason',string(reason),'inspection',inspection);
end

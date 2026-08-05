function receipt = run_Dual_analysis3_rec(request, runtime)
%RUN_DUAL_ANALYSIS3_REC Schedule Dual_analysis3 once per selected Cycle.
% The Record layer discovers jobs, resolves shared ROI/reuse sources, records
% status, and optionally invokes Dual_rec_analysis3. All single-Cycle
% scientific processing is delegated to aaa.workflows.Dual_analysis3.

if nargin < 1 || isempty(request)
    request = aaa.helpers.common.direct_run_request("dual", "record");
end
if nargin < 2 || isempty(runtime)
    runtime = struct();
end

request = apply_record_checkbox_overlays(request);
[request, validation] = aaa.schema.validate_request(request);
if ~validation.valid
    error('AAA:DualRecord:InvalidRequest', '%s', ...
        strjoin(string(validation.errors), newline));
end
execution_plan = validation.execution_plan;
[~,preflight_match]=aaa.helpers.common.runtime_preflight(runtime,request);
if ~request.record.record_average_only&&preflight_match~="exact"
    input_inspection=aaa.helpers.common.inspect_input_structure(request,execution_plan);
    if ~input_inspection.valid
        error('AAA:DualRecord:InvalidInputStructure','%s', ...
            strjoin(input_inspection.errors,newline));
    end
end
rec_path = char(string(request.input_path));

receipt = empty_record_receipt(request);
[record_runtime, record_logger, record_log_file, owns_record_logger] = ...
    aaa.io.attach_analysis_logger(runtime, request, rec_path);
record_runtime.aaa_cycle_runtime = runtime;
emit_record_event(record_runtime, request, "", "workflow", "running", ...
    "Starting Dual Record workflow.", rec_path);
if request.workflow.write_manifest
    [~, manifest_file] = aaa.io.create_manifest(request, execution_plan, rec_path);
    receipt.manifest_file = string(manifest_file);
    aaa.io.update_manifest(receipt.manifest_file, struct( ...
        'status', "running", 'log_file', string(record_log_file), ...
        'started_at', datetime('now')));
end

diary_cleanup = [];
if request.record.write_diary
    diary_file = fullfile(rec_path, sprintf('Dual_analysis3_rec_%s.log', ...
        char(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss'))));
    diary(diary_file);
    diary_cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
    fprintf('[AAA] Dual Record diary: %s\n', diary_file);
end

try
    cycle_entries = aaa.helpers.common.discover_record_cycles( ...
        rec_path, request.record.cycles, true);
    batch_summary_file = fullfile(rec_path, 'Dual_analysis3_rec_batch_summary.mat');
    analysis_mode = char(request.workflow.preset); %#ok<NASGU>
    reference_roi_file = "";
    job_plan = repmat(empty_job(), 0, 1);

    if request.record.record_average_only
        [batch_results, job_plan, reference_roi_file] = ...
            scan_existing_record_results(cycle_entries, request);
    elseif request.workflow.preset == "motion_only"
        [batch_results, job_plan] = run_motion_jobs( ...
            cycle_entries, request, record_runtime);
    elseif request.workflow.preset == "analysis_only"
        [batch_results, job_plan] = run_analysis_only_jobs( ...
            cycle_entries, request, execution_plan, record_runtime);
        reference_roi_file = first_completed_roi(batch_results);
    else
        [batch_results, job_plan, reference_roi_file] = run_shared_roi_jobs( ...
            cycle_entries, request, execution_plan, record_runtime);
    end

    analysis_only_output_root_dir = char(string(request.record.output_root)); %#ok<NASGU>
    analysis_only_run_tag = string(request.record.run_tag); %#ok<NASGU>
    voltage_polarity = request.params.dual.voltage_polarity; %#ok<NASGU>
    calcium_polarity = request.params.dual.calcium_polarity; %#ok<NASGU>
    calcium_smoothing_window = request.params.dual.calcium_smoothing_window; %#ok<NASGU>
    dual3_rec_control_effective = request; %#ok<NASGU>
    save(batch_summary_file, 'batch_results', 'reference_roi_file', ...
        'analysis_mode', 'analysis_only_output_root_dir', 'analysis_only_run_tag', ...
        'voltage_polarity', 'calcium_polarity', 'calcium_smoothing_window', ...
        'dual3_rec_control_effective', 'job_plan');

    record_average_receipt = struct();
    if request.record.record_average
        average_control = struct( ...
            'mode', "dual", ...
            'record_path', string(rec_path), ...
            'cycle_results', batch_results, ...
            'reference_roi_file', string(reference_roi_file), ...
            'source_results_path', string(request.workflow.source_results_path), ...
            'output_dir_name', "Dual_analysis3_rec_average", ...
            'request', request, ...
            'voltage_polarity', request.params.dual.voltage_polarity, ...
            'calcium_polarity', request.params.dual.calcium_polarity, ...
            'calcium_smoothing_window', request.params.dual.calcium_smoothing_window);
        record_average_receipt = Record_analysis3(average_control, record_runtime);
    end

    receipt.status = record_completion_status(batch_results);
    receipt.output_path = string(rec_path);
    receipt.batch_summary_file = string(batch_summary_file);
    receipt.batch_results = batch_results;
    receipt.job_plan = job_plan;
    receipt.reference_roi_file = string(reference_roi_file);
    receipt.record_average = record_average_receipt;
    receipt.completed_sections = completed_section_names(execution_plan);
    receipt.output_files = string(batch_summary_file);
    if isstruct(record_average_receipt) && isfield(record_average_receipt, 'output_files')
        receipt.output_files = unique([receipt.output_files; ...
            string(record_average_receipt.output_files(:))], 'stable');
    end
    if strlength(receipt.manifest_file) > 0
        aaa.io.update_manifest(receipt.manifest_file, struct( ...
            'status', receipt.status, ...
            'completed_sections', receipt.completed_sections, ...
            'failed_section', receipt.failed_section, ...
            'output_files', receipt.output_files, ...
            'log_file', string(record_log_file), ...
            'completed_at', datetime('now'), ...
            'error', receipt.error, ...
            'record', struct('batch_results', batch_results, 'job_plan', job_plan, ...
                'reference_roi_file', string(reference_roi_file))));
    end
    emit_record_event(record_runtime, request, "", "workflow", receipt.status, ...
        "Dual Record workflow finished.", rec_path);
    close_owned_logger(record_logger, owns_record_logger);
catch ME
    receipt.status = "failed";
    receipt.failed_section = "record_scheduler";
    receipt.error = exception_struct(ME);
    if strlength(receipt.manifest_file) > 0
        aaa.io.update_manifest(receipt.manifest_file, struct( ...
            'status', receipt.status, ...
            'failed_section', receipt.failed_section, ...
            'log_file', string(record_log_file), ...
            'completed_at', datetime('now'), ...
            'error', receipt.error));
    end
    emit_record_event(record_runtime, request, "", "workflow", "failed", ...
        string(ME.message), rec_path);
    close_owned_logger(record_logger, owns_record_logger);
    rethrow(ME);
end
end

function request = apply_record_checkbox_overlays(request)
% Preserve the old Record checkboxes as explicit UI controls while keeping
% the section table authoritative when the user set an overlay directly.
if ~isfield(request, 'workflow') || ~isfield(request.workflow, 'sections') ...
        || isempty(request.workflow.sections)
    request.workflow.sections = struct();
end
if isfield(request, 'record') && isfield(request.record, 'run_motion') ...
        && isfield(request, 'workflow') && isfield(request.workflow, 'preset') ...
        && strcmpi(string(request.workflow.preset), "full") ...
        && ~isfield(request.workflow.sections, 'motion')
    request.workflow.sections.motion = ...
        choose(logical(request.record.run_motion), "run", "skip");
end
if isfield(request, 'record') && isfield(request.record, 'reuse_saved_peaks') ...
        && isfield(request.workflow, 'preset') ...
        && strcmpi(string(request.workflow.preset), "analysis_only") ...
        && ~isfield(request.workflow.sections, 'peak')
    request.workflow.sections.peak = ...
        choose(logical(request.record.reuse_saved_peaks), "reuse", "run");
end
end

function [batch_results, job_plan] = run_motion_jobs(entries, request, runtime)
batch_results = repmat(empty_batch_result(), 0, 1);
job_plan = repmat(empty_job(), 0, 1);
for idx = 1:numel(entries)
    entry = entries(idx);
    if request.record.skip_existing
        existing = find_latest_motion_result(entry.cycle_path);
        if strlength(existing) > 0
            job = make_job(entry, "skip_existing", "", "", existing);
            job_plan(end+1,1) = job; %#ok<AGROW>
            batch_results(end+1,1) = make_batch_result( ...
                entry, "skipped_existing_motion_correction", ...
                string(fileparts(existing)), "", existing); %#ok<AGROW>
            emit_record_event(runtime, request, entry.label, "cycle", "skipped", ...
                "Existing motion output reused.", fileparts(existing));
            continue;
        end
    end
    cycle_request = prepare_cycle_request(request, entry);
    cycle_request.workflow.preset = "motion_only";
    cycle_request.workflow.sections = request.workflow.sections;
    job = make_job(entry, "run", cycle_request.workflow.preset, "", "");
    job_plan(end+1,1) = job; %#ok<AGROW>
    [result, ok] = execute_cycle(cycle_request, entry, request, runtime);
    batch_results(end+1,1) = result; %#ok<AGROW>
    if ~ok && request.record.stop_on_error
        error('AAA:DualRecord:CycleFailed', '%s failed: %s', entry.label, result.message);
    end
end
end

function [batch_results, job_plan] = run_analysis_only_jobs(entries, request, plan, runtime)
batch_results = repmat(empty_batch_result(), 0, 1);
job_plan = repmat(empty_job(), 0, 1);
run_tag = string(request.record.run_tag);
if strlength(run_tag) == 0
    run_tag = string(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss'));
end
for idx = 1:numel(entries)
    entry = entries(idx);
    probe=request;probe.scope="single";probe.input_path=string(entry.cycle_path);
    [~,match_kind,preflight_cycle]=aaa.helpers.common.runtime_preflight(runtime,probe);
    source_path="";source_report=struct('status',"cached");
    if match_kind=="cycle"
        source_path=string(preflight_cycle.selected_reuse_source);
    end
    if strlength(source_path)==0
        search_root = string(request.workflow.source_results_path);
        if strlength(search_root) == 0, search_root = string(entry.cycle_path); end
        [source_path, source_report] = aaa.helpers.common.resolve_reuse_source( ...
            "dual", search_root, plan, ...
            struct('cycle_name',entry.cycle_name,'label',entry.label));
    end
    if strlength(source_path) == 0
        batch_results(end+1,1) = make_batch_result( ...
            entry, "missing_source_result", "", "", ...
            "No existing Dual_analysis3 result satisfies analysis_only reuse (" + ...
            source_report.status + ")."); %#ok<AGROW>
        job_plan(end+1,1) = make_job(entry, "missing_source", "analysis_only", "", ""); %#ok<AGROW>
        continue;
    end
    if request.record.skip_existing
        existing = find_latest_analysis_only_result(entry.cycle_path, request.record.output_root);
        if strlength(existing) > 0
            batch_results(end+1,1) = make_batch_result( ...
                entry, "skipped_existing_analysis_only_result", fileparts(existing), ...
                resolve_roi_path_in_result(fileparts(existing)), existing); %#ok<AGROW>
            job_plan(end+1,1) = make_job(entry, "skip_existing", "analysis_only", source_path, existing); %#ok<AGROW>
            continue;
        end
    end
    cycle_request = prepare_cycle_request(request, entry);
    cycle_request.workflow.preset = "analysis_only";
    cycle_request.workflow.sections = request.workflow.sections;
    cycle_request.workflow.source_results_path = source_path;
    cycle_request.workflow.run_name = make_cycle_run_name(entry, run_tag, request.workflow.run_name);
    cycle_request.workflow.output_path = fullfile(entry.cycle_path, ...
        char(string(request.record.output_root)), char(cycle_request.workflow.run_name));
    job_plan(end+1,1) = make_job(entry, "run", "analysis_only", source_path, ""); %#ok<AGROW>
    [result, ok] = execute_cycle(cycle_request, entry, request, runtime);
    batch_results(end+1,1) = result; %#ok<AGROW>
    if ~ok && request.record.stop_on_error
        error('AAA:DualRecord:CycleFailed', '%s failed: %s', entry.label, result.message);
    end
end
end

function [batch_results, job_plan, reference_roi_file] = ...
        run_shared_roi_jobs(entries, request, plan, runtime)
batch_results = repmat(empty_batch_result(), 0, 1);
job_plan = repmat(empty_job(), 0, 1);
reference_roi_file = resolve_explicit_reference_roi(request);
reference_idx = find(strcmpi(string({entries.cycle_name}), ...
    string(request.record.reference_cycle_name)), 1, 'first');
if isempty(reference_idx)
    error('AAA:DualRecord:MissingReferenceCycle', ...
        'Reference Cycle was not found: %s', request.record.reference_cycle_name);
end
reference_entry = entries(reference_idx);

if strlength(reference_roi_file) == 0 ...
        && request.record.reuse_reference_roi ...
        && ~request.record.redraw_reference_roi
    reference_roi_file = find_latest_roi_file(reference_entry.cycle_path);
end

reference_was_run = false;
if strlength(reference_roi_file) == 0 || request.record.redraw_reference_roi
    if plan.roi ~= "run"
        error('AAA:DualRecord:ReferenceRoiUnavailable', ...
            ['No reusable reference ROI was found, while the resolved plan has roi="%s". ' ...
             'Choose an ROI path or a plan that runs ROI.'], plan.roi);
    end
    if request.record.reference_roi_only && request.record.include_reference_cycle
        error('AAA:DualRecord:ReferencePolicyConflict', ...
            ['reference_roi_only and include_reference_cycle cannot both be selected: ' ...
             'each Cycle may call Dual_analysis3 at most once.']);
    end
    ref_request = prepare_cycle_request(request, reference_entry);
    ref_request.workflow.offset_mode = request.record.reference_offset_mode;
    ref_request.workflow.roi_path = "";
    ref_request = resolve_cycle_source(ref_request, plan, "", reference_entry,runtime);
    if request.record.reference_roi_only
        ref_request.workflow.preset = "custom";
        ref_request.workflow.sections = roi_only_sections(plan);
    end
    job_plan(end+1,1) = make_job(reference_entry, ...
        choose(request.record.reference_roi_only, "reference_roi_only", "reference_full"), ...
        ref_request.workflow.preset, ref_request.workflow.source_results_path, ""); %#ok<AGROW>
    [ref_result, ok] = execute_cycle(ref_request, reference_entry, request, runtime);
    batch_results(end+1,1) = ref_result; %#ok<AGROW>
    if ~ok
        error('AAA:DualRecord:ReferenceCycleFailed', ...
            'Reference Cycle failed and the shared ROI is unavailable: %s', ref_result.message);
    end
    reference_roi_file = string(ref_result.roi_file);
    if strlength(reference_roi_file) == 0 || ~isfile(reference_roi_file)
        error('AAA:DualRecord:ReferenceRoiMissing', ...
            'Reference Cycle completed without a valid ROI output: %s', reference_entry.cycle_path);
    end
    reference_was_run = true;
end

roi_folder = string(fileparts(char(reference_roi_file)));
for idx = 1:numel(entries)
    entry = entries(idx);
    is_reference = idx == reference_idx;
    if is_reference && reference_was_run
        continue;
    end
    if is_reference && ~request.record.include_reference_cycle
        batch_results(end+1,1) = make_batch_result( ...
            entry, "skipped_reference_already_used_for_roi", "", ...
            reference_roi_file, "Reference Cycle supplied the shared ROI."); %#ok<AGROW>
        job_plan(end+1,1) = make_job(entry, "skip_reference", request.workflow.preset, "", reference_roi_file); %#ok<AGROW>
        continue;
    end
    if request.record.skip_existing
        existing = find_latest_explicit_result(entry.cycle_path);
        if strlength(existing) > 0
            batch_results(end+1,1) = make_batch_result( ...
                entry, "skipped_existing_result", fileparts(existing), ...
                reference_roi_file, existing); %#ok<AGROW>
            job_plan(end+1,1) = make_job(entry, "skip_existing", request.workflow.preset, "", existing); %#ok<AGROW>
            continue;
        end
    end
    cycle_request = prepare_cycle_request(request, entry);
    cycle_request.workflow.sections = request.workflow.sections;
    cycle_request.workflow.sections.roi = "reuse";
    cycle_request.workflow.sections.registration = "skip";
    cycle_request.workflow.roi_path = roi_folder;
    cycle_request.workflow.offset_mode = request.record.cycle_offset_mode;
    [~, ~, ~, cycle_plan] = aaa.schema.resolve_sections( ...
        "dual", cycle_request.workflow.preset, cycle_request.workflow.sections);
    cycle_request = resolve_cycle_source(cycle_request, cycle_plan, roi_folder, entry,runtime);
    job_plan(end+1,1) = make_job(entry, "run", cycle_request.workflow.preset, ...
        cycle_request.workflow.source_results_path, reference_roi_file); %#ok<AGROW>
    [result, ok] = execute_cycle(cycle_request, entry, request, runtime);
    batch_results(end+1,1) = result; %#ok<AGROW>
    if ~ok && request.record.stop_on_error
        error('AAA:DualRecord:CycleFailed', '%s failed: %s', entry.label, result.message);
    end
end
end

function [batch_results, job_plan, reference_roi_file] = ...
        scan_existing_record_results(entries, request)
batch_results = repmat(empty_batch_result(), 0, 1);
job_plan = repmat(empty_job(), 0, 1);
for idx = 1:numel(entries)
    entry = entries(idx);
    existing = find_latest_explicit_result(entry.cycle_path);
    if strlength(existing) > 0
        result_dir = string(fileparts(existing));
        batch_results(end+1,1) = make_batch_result( ...
            entry, "existing_result_detected", result_dir, ...
            resolve_roi_path_in_result(result_dir), existing); %#ok<AGROW>
        job_plan(end+1,1) = make_job(entry, "existing_only", ...
            request.workflow.preset, result_dir, existing); %#ok<AGROW>
    else
        batch_results(end+1,1) = make_batch_result( ...
            entry, "missing_result", "", "", ...
            "No completed Dual result with a reusable trace was found under this Cycle."); %#ok<AGROW>
        job_plan(end+1,1) = make_job(entry, "missing_result", ...
            request.workflow.preset, "", ""); %#ok<AGROW>
    end
end
reference_match = find(strcmpi(string({entries.cycle_name}), ...
    string(request.record.reference_cycle_name)), 1, 'first');
if isempty(reference_match)
    reference_roi_file = "";
else
    reference_roi_file = find_latest_roi_file(entries(reference_match).cycle_path);
end
end

function cycle_request = prepare_cycle_request(record_request, entry)
cycle_request = record_request;
cycle_request.scope = "single";
cycle_request.input_path = string(entry.cycle_path);
cycle_request.record = struct();
cycle_request.workflow.output_path = "";
cycle_request.workflow.run_name = make_cycle_run_name( ...
    entry, string(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss-SSS')), ...
    record_request.workflow.run_name);
end

function cycle_request = resolve_cycle_source(cycle_request, plan, roi_path, entry,runtime)
if ~plan_requires_source(plan, roi_path)
    return;
end
search_root = string(cycle_request.workflow.source_results_path);
if strlength(search_root) == 0, search_root = string(cycle_request.input_path); end
[~,match_kind,preflight_cycle]=aaa.helpers.common.runtime_preflight(runtime,cycle_request);
source_path="";report=struct('status',"cached");
if match_kind=="cycle"
    source_path=string(preflight_cycle.selected_reuse_source);
end
if strlength(source_path)==0
    [source_path, report] = aaa.helpers.common.resolve_reuse_source( ...
        "dual", search_root, plan, ...
        struct('cycle_name',entry.cycle_name,'label',entry.label));
end
if strlength(source_path) == 0
    error('AAA:DualRecord:ReuseSourceMissing', ...
        'No result folder satisfies the reuse plan for %s under: %s (%s)', ...
        entry.label, search_root, report.status);
end
cycle_request.workflow.source_results_path = source_path;
end

function [batch_result, ok] = execute_cycle(cycle_request, entry, record_request, runtime)
emit_record_event(runtime, record_request, entry.label, "cycle", "running", ...
    "Calling Dual_analysis3 once; reuse_source=" + ...
    string(cycle_request.workflow.source_results_path), entry.cycle_path);
try
    cycle_runtime = runtime;
    if isfield(runtime, 'aaa_cycle_runtime')
        cycle_runtime = runtime.aaa_cycle_runtime;
    end
    cycle_receipt = aaa.workflows.Dual_analysis3(cycle_request, cycle_runtime);
    batch_result = make_batch_result(entry, "completed", ...
        cycle_receipt.output_path, cycle_receipt.roi_file, "");
    ok = true;
    emit_record_event(runtime, record_request, entry.label, "cycle", "completed", ...
        "Cycle completed.", cycle_receipt.output_path);
catch ME
    report = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
    error_log = save_cycle_error_report(entry.record_path, entry.label, report);
    batch_result = make_batch_result(entry, "failed", "", "", ...
        report + newline + "Error log: " + error_log);
    ok = false;
    emit_record_event(runtime, record_request, entry.label, "cycle", "failed", ...
        string(ME.message), entry.cycle_path);
end
end

function sections = roi_only_sections(plan)
names = fieldnames(plan);
sections = struct();
for idx = 1:numel(names)
    sections.(names{idx}) = "skip";
end
sections.input = "run";
sections.motion = plan.motion;
sections.registration = plan.registration;
sections.roi = "run";
end

function roi_file = resolve_explicit_reference_roi(request)
path_value = string(request.record.reference_roi_path);
if strlength(path_value) == 0
    path_value = string(request.workflow.roi_path);
end
if strlength(path_value) == 0
    roi_file = "";
elseif isfile(path_value)
    assert_valid_roi_file(path_value);
    roi_file = path_value;
elseif isfolder(path_value)
    roi_file = resolve_roi_file_in_folder(path_value);
else
    error('AAA:DualRecord:InvalidReferenceRoiPath', ...
        'Reference ROI path does not exist: %s', path_value);
end
end

function roi_file = resolve_roi_file_in_folder(folder_path)
priority = [fullfile("results","roi.mat"), ...
    "1_dual_roi_results.mat","1_raw_ROI.mat","roi.mat"];
for idx = 1:numel(priority)
    candidate = fullfile(char(folder_path), char(priority(idx)));
    if isfile(candidate) && is_valid_roi_file(candidate)
        roi_file = string(candidate);
        return;
    end
end
listing = dir(fullfile(char(folder_path), '*ROI*.mat'));
valid = strings(0,1);
for idx = 1:numel(listing)
    candidate = string(fullfile(listing(idx).folder, listing(idx).name));
    if is_valid_roi_file(candidate)
        valid(end+1,1) = candidate; %#ok<AGROW>
    end
end
if isempty(valid)
    error('AAA:DualRecord:NoValidRoi', ...
        'No valid ROI MAT file was found directly in: %s', folder_path);
elseif numel(valid) > 1
    error('AAA:DualRecord:AmbiguousRoi', ...
        'Multiple valid ROI MAT files were found in: %s', folder_path);
end
roi_file = valid(1);
end

function roi_file = find_latest_roi_file(cycle_path)
patterns = {fullfile('results','roi.mat'),'1_dual_roi_results.mat', ...
    '1_raw_ROI.mat','roi.mat','*ROI*.mat'};
listing = struct([]);
for idx = 1:numel(patterns)
    found = dir(fullfile(char(cycle_path), '**', patterns{idx}));
    listing = [listing; found(~[found.isdir])]; %#ok<AGROW>
end
roi_file = "";
if isempty(listing), return; end
[~, order] = sort([listing.datenum], 'descend');
for idx = order(:)'
    candidate = string(fullfile(listing(idx).folder, listing(idx).name));
    if is_valid_roi_file(candidate)
        roi_file = candidate;
        return;
    end
end
end

function assert_valid_roi_file(file_path)
if ~is_valid_roi_file(file_path)
    error('AAA:DualRecord:InvalidRoiFile', 'Invalid ROI MAT file: %s', file_path);
end
end

function tf = is_valid_roi_file(file_path)
tf = false;
try
    [results,~]=load_analysis_results(file_path,struct('mode',"dual"));
    if isfield(results,'roi')&&isstruct(results.roi)
        value=results.roi;
        if isfield(value,'rois'),value=value.rois;end
        tf=(isfield(value,'bwmask')&&~isempty(value.bwmask)) ...
            ||(isfield(value,'voltage')&&~isempty(value.voltage));
        if tf,return;end
    end
catch
end
try,saved=load(char(file_path));catch,return;end
tf=(isfield(saved,'rois')&&isstruct(saved.rois) ...
        &&isfield(saved.rois,'bwmask')&&~isempty(saved.rois.bwmask)) ...
    ||(isfield(saved,'bwmask')&&~isempty(saved.bwmask)) ...
    ||(isfield(saved,'mask')&&~isempty(saved.mask));
end

function marker = find_latest_explicit_result(cycle_path)
marker = find_latest_completed_result(cycle_path,"Dual_analysis3","");
end

function marker = find_latest_motion_result(cycle_path)
listing = dir(fullfile(char(cycle_path), 'Dual_analysis3', ...
    '**', 'shared_motion_shifts_result.mat'));
[~, order] = sort([listing.datenum], 'descend');
marker = "";
for idx = order(:)'
    folder = listing(idx).folder;
    if aaa.io.manifest_allows_reuse(folder) ...
            && ~isempty(dir(fullfile(folder,'voltage_motion_corrected_ds*.tif'))) ...
            && ~isempty(dir(fullfile(folder,'calcium_motion_corrected_ds*.tif')))
        marker = string(fullfile(folder, listing(idx).name));
        return;
    end
end
end

function marker = find_latest_analysis_only_result(cycle_path, output_root)
marker=find_latest_completed_result(cycle_path,output_root,"analysis_only");
end

function marker=find_latest_completed_result(cycle_path,output_root,preset)
% Current AAA results are identified by a completed manifest plus a trace
% result record. Historical explicit marker files remain read-only fallback
% candidates and participate in the same newest-first selection.
root=fullfile(char(cycle_path),char(string(output_root)));
paths=strings(0,1);timestamps=zeros(0,1);
manifests=dir(fullfile(root,'**','analysis_manifest.mat'));
for idx=1:numel(manifests)
    file=string(fullfile(manifests(idx).folder,manifests(idx).name));
    try
        manifest=aaa.io.load_manifest(file);
        if ~manifest_is_completed_trace_result(manifest,preset),continue;end
        paths(end+1,1)=file; %#ok<AGROW>
        timestamps(end+1,1)=manifests(idx).datenum; %#ok<AGROW>
    catch
    end
end
legacy=dir(fullfile(root,'**','-1_explicit_dual_results.mat'));
for idx=1:numel(legacy)
    if strlength(preset)>0&&~legacy_result_matches_preset(legacy(idx).folder,preset)
        continue;
    end
    if ~aaa.io.manifest_allows_reuse(legacy(idx).folder),continue;end
    paths(end+1,1)=string(fullfile(legacy(idx).folder,legacy(idx).name)); %#ok<AGROW>
    timestamps(end+1,1)=legacy(idx).datenum; %#ok<AGROW>
end
marker="";
if isempty(paths),return;end
[~,order]=sort(timestamps,'descend');
marker=paths(order(1));
end

function tf=manifest_is_completed_trace_result(manifest,preset)
tf=isfield(manifest,'status')&&strcmpi(string(manifest.status),"completed");
if ~tf,return;end
if strlength(preset)>0
    tf=isfield(manifest,'request')&&isstruct(manifest.request) ...
        &&isfield(manifest.request,'workflow') ...
        &&isfield(manifest.request.workflow,'preset') ...
        &&strcmpi(string(manifest.request.workflow.preset),preset);
    if ~tf,return;end
end
tf=false;
if ~isfield(manifest,'results')||~isstruct(manifest.results),return;end
for idx=1:numel(manifest.results)
    item=manifest.results(idx);
    if isfield(item,'section')&&strcmpi(string(item.section),"trace") ...
            &&isfield(item,'status') ...
            &&ismember(lower(string(item.status)),["completed","reused"])
        file="";
        if isfield(item,'file'),file=string(item.file);end
        if strlength(file)>0&&~isfile(file)&&isfield(manifest,'output_path')
            relative="";
            if isfield(item,'relative_file'),relative=string(item.relative_file);end
            if strlength(relative)>0,file=fullfile(string(manifest.output_path),relative);end
        end
        tf=strlength(file)>0&&isfile(file);
        return;
    end
end
end

function tf=legacy_result_matches_preset(folder,preset)
tf=false;
try
    saved=load(fullfile(folder,'dual_info.mat'),'dual_info');
    tf=isfield(saved,'dual_info') ...
        &&isfield(saved.dual_info,'workflow_control') ...
        &&isfield(saved.dual_info.workflow_control,'preset') ...
        &&strcmpi(string(saved.dual_info.workflow_control.preset),preset);
catch
end
end

function tf = plan_requires_source(plan, roi_path)
source_plan = plan;
if source_plan.roi == "reuse" && strlength(string(roi_path)) > 0
    source_plan.roi = "skip";
end
tf = any(structfun(@(value) string(value) == "reuse", source_plan));
end

function value = latest_file(listing)
value = "";
if isempty(listing), return; end
[~, order] = sort([listing.datenum], 'descend');
for idx = order(:)'
    if aaa.io.manifest_allows_reuse(listing(idx).folder)
        value = string(fullfile(listing(idx).folder, listing(idx).name));
        return;
    end
end
end

function roi = resolve_roi_path_in_result(result_dir)
candidates=[string(fullfile(result_dir,'results','roi.mat')); ...
    string(fullfile(result_dir,'1_dual_roi_results.mat'))];
roi="";for idx=1:numel(candidates),if isfile(candidates(idx)),roi=candidates(idx);return;end,end
end

function name = make_cycle_run_name(entry, tag, explicit_name)
if strlength(string(explicit_name)) > 0
    name = string(explicit_name);
else
    name = entry.record_name + "_" + entry.cycle_name + "_" + string(tag);
end
name = regexprep(name, '[<>:"/\\|?*]', '_');
end

function result = make_batch_result(entry, status, save_path, roi_file, message)
result = empty_batch_result();
result.cycle_name = string(entry.label);
result.cycle_path = string(entry.cycle_path);
result.status = string(status);
result.save_path = string(save_path);
result.roi_file = string(roi_file);
result.message = string(message);
end

function result = empty_batch_result()
result = struct('cycle_name',"",'cycle_path',"",'status',"", ...
    'save_path',"",'roi_file',"",'message',"");
end

function job = make_job(entry, action, preset, source_path, context_path)
job = empty_job();
job.cycle_name = string(entry.label);
job.cycle_path = string(entry.cycle_path);
job.action = string(action);
job.preset = string(preset);
job.source_results_path = string(source_path);
job.context_path = string(context_path);
end

function job = empty_job()
job = struct('cycle_name',"",'cycle_path',"",'action',"", ...
    'preset',"",'source_results_path',"",'context_path',"");
end

function receipt = empty_record_receipt(request)
receipt = struct('status',"pending",'mode',"dual",'scope',"record", ...
    'input_path',string(request.input_path),'output_path',"", ...
    'manifest_file',"",'completed_sections',strings(0,1), ...
    'failed_section',"",'error',struct(),'output_files',strings(0,1), ...
    'batch_summary_file',"",'batch_results',repmat(empty_batch_result(),0,1), ...
    'job_plan',repmat(empty_job(),0,1),'reference_roi_file',"", ...
    'record_average',struct());
end

function status = record_completion_status(results)
if any(strcmpi(string({results.status}), "failed"))
    status = "completed_with_errors";
else
    status = "completed";
end
end

function roi_file = first_completed_roi(results)
roi_file = "";
for idx = 1:numel(results)
    if strlength(string(results(idx).roi_file)) > 0
        roi_file = string(results(idx).roi_file);
        return;
    end
end
end

function names = completed_section_names(plan)
all_names = string(fieldnames(plan));
keep = false(size(all_names));
for idx = 1:numel(all_names)
    keep(idx) = string(plan.(all_names(idx))) ~= "skip";
end
names = all_names(keep);
end

function log_file = save_cycle_error_report(record_path, cycle_name, report)
error_dir = fullfile(char(record_path), 'Dual_analysis3_rec_error_logs');
if ~isfolder(error_dir), mkdir(error_dir); end
safe_name = regexprep(string(cycle_name), '[<>:"/\\|?*]', '_');
log_file = string(fullfile(error_dir, sprintf('%s_%s_error.txt', ...
    char(safe_name), char(datetime('now','Format','yyyy-MM-dd HH-mm-ss-SSS')))));
fid = fopen(log_file, 'w');
if fid < 0
    warning('AAA:DualRecord:ErrorLogWriteFailed', ...
        'Could not write Cycle error log: %s', log_file);
    return;
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', report);
end

function emit_record_event(runtime, request, cycle, section, state, message, output_path)
event = struct('timestamp',datetime('now'),'type',"section", ...
    'mode',"dual",'scope',"record",'input_path',string(request.input_path), ...
    'cycle',string(cycle),'section',string(section),'state',string(state), ...
    'message',string(message),'output_path',string(output_path),'details',struct());
aaa.io.emit_event(runtime, event);
end

function close_owned_logger(logger, owns_logger)
if owns_logger && ~isempty(logger) && isa(logger, 'aaa.io.AnalysisLogger')
    logger.close();
end
end

function value = choose(condition, true_value, false_value)
if condition, value = true_value; else, value = false_value; end
end

function result = exception_struct(ME)
result = struct('identifier',string(ME.identifier),'message',string(ME.message), ...
    'report',string(getReport(ME,'extended','hyperlinks','off')));
end

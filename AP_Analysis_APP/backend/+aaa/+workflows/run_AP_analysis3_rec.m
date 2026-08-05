function receipt = run_AP_analysis3_rec(request, runtime)
%RUN_AP_ANALYSIS3_REC Schedule AP_analysis3 once for each selected Cycle.
% This first AP Record backend performs Cycle discovery, optional shared-ROI
% scheduling, error handling, logging, and status persistence only. It does
% not calculate a Record average.

if nargin < 1 || isempty(request)
    request = aaa.helpers.common.direct_run_request("ap", "record");
end
if nargin < 2 || isempty(runtime)
    runtime = struct();
end

[request, validation] = aaa.schema.validate_request(request);
if ~validation.valid
    error('AAA:APRecord:InvalidRequest', '%s', strjoin(string(validation.errors), newline));
end
if request.mode ~= "ap" || request.scope ~= "record"
    error('AAA:APRecord:WrongRequest', ...
        'run_AP_analysis3_rec requires mode="ap" and scope="record".');
end
 [~,preflight_match]=aaa.helpers.common.runtime_preflight(runtime,request);
if ~request.record.record_average_only&&preflight_match~="exact"
    input_inspection=aaa.helpers.common.inspect_input_structure( ...
        request,validation.execution_plan);
    if ~input_inspection.valid
        error('AAA:APRecord:InvalidInputStructure','%s', ...
            strjoin(input_inspection.errors,newline));
    end
end
if request.record.reference_roi_only
    error('AAA:APRecord:ReferenceROIOnlyNotMigrated', ...
        ['reference_roi_only requires a stage-addressable AP pipeline. The ' ...
         'current frozen AP bridge can only run the complete downstream workflow.']);
end

execution_plan = validation.execution_plan;
items = aaa.helpers.common.discover_record_cycles( ...
    request.input_path, request.record.cycles, true);
if isempty(items)
    error('AAA:APRecord:NoCycles', ...
        'No selected Cycle folders were found under %s.', request.input_path);
end

run_tag = string(request.record.run_tag);
if strlength(run_tag) == 0
    run_tag = string(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss'));
end
summary_path = fullfile(request.input_path, 'AP_analysis3_rec', char(run_tag));
if ~isfolder(summary_path)
    mkdir(summary_path);
end
[runtime, analysis_logger, log_file, owns_logger] = ...
    aaa.io.attach_analysis_logger(runtime, request, summary_path);
logger_cleanup = onCleanup(@() close_owned_logger(analysis_logger, owns_logger)); %#ok<NASGU>
[~, manifest_file] = aaa.io.create_manifest( ...
    request, execution_plan, summary_path);
update_manifest_if_present(manifest_file, struct('log_file', string(log_file)));

receipt = empty_record_receipt(request, summary_path, manifest_file);
emit_event_safe(runtime, make_event( ...
    "record_started", request, summary_path, ...
    sprintf('Scheduling %d AP Cycle job(s).', numel(items))));

diary_cleanup = [];
diary_file = "";
if request.record.write_diary
    diary_file = string(fullfile(summary_path, 'AP_analysis3_rec_diary.txt'));
    diary(char(diary_file));
    diary on;
    diary_cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
end

cycle_receipts = repmat(empty_cycle_receipt(), 0, 1);
reference_roi_file = "";
if ~request.record.record_average_only
    reference_roi_file = resolve_initial_reference_roi(request, items, execution_plan);
end
reference_cycle_idx = find(string({items.cycle_name}) == ...
    string(request.record.reference_cycle_name), 1, 'first');
reference_has_run = false;

try
    if ~request.record.record_average_only
    needs_shared_roi = execution_plan.roi == "run" && ( ...
        request.record.reuse_reference_roi ...
        || request.record.redraw_reference_roi ...
        || strlength(string(request.record.reference_roi_path)) > 0);

    if needs_shared_roi && strlength(reference_roi_file) == 0
        if isempty(reference_cycle_idx)
            error('AAA:APRecord:ReferenceCycleNotFound', ...
                'Reference Cycle was not found: %s', request.record.reference_cycle_name);
        end
        reference_request = build_cycle_request( ...
            request, items(reference_cycle_idx), run_tag, "",runtime);
        reference_request.workflow.sections.roi = "run";
        reference_request.workflow.roi_path = "";
        reference_receipt = run_one_cycle(reference_request, runtime, items(reference_cycle_idx));
        cycle_receipts(end + 1, 1) = reference_receipt; %#ok<AGROW>
        reference_has_run = true;
        reference_roi_file = string(reference_receipt.roi_file);
        if strlength(reference_roi_file) == 0 || ~isfile(reference_roi_file)
            error('AAA:APRecord:ReferenceROINotProduced', ...
                'Reference Cycle completed without producing results/roi.mat.');
        end
    end

    for item_idx = 1:numel(items)
        is_reference = item_idx == reference_cycle_idx;
        if is_reference && reference_has_run
            continue;
        end
        if is_reference && ~request.record.include_reference_cycle ...
                && strlength(reference_roi_file) > 0
            continue;
        end

        cycle_roi_file = "";
        if execution_plan.roi == "reuse"
            cycle_roi_file = string(request.workflow.roi_path);
            if strlength(cycle_roi_file) == 0
                cycle_roi_file = string(request.workflow.source_results_path);
            end
        elseif needs_shared_roi && strlength(reference_roi_file) > 0
            cycle_roi_file = reference_roi_file;
        end
        cycle_request = build_cycle_request(request, items(item_idx), run_tag, cycle_roi_file,runtime);

        existing_receipt = maybe_existing_receipt(cycle_request, items(item_idx));
        if request.record.skip_existing && existing_receipt.status == "skipped"
            cycle_receipts(end + 1, 1) = existing_receipt; %#ok<AGROW>
            continue;
        end

        try
            cycle_receipts(end + 1, 1) = ...
                run_one_cycle(cycle_request, runtime, items(item_idx)); %#ok<AGROW>
        catch ME
            cycle_receipts(end + 1, 1) = ...
                failed_cycle_receipt(items(item_idx), cycle_request, ME); %#ok<AGROW>
            if request.record.stop_on_error
                rethrow(ME);
            end
        end
    end
    end

    receipt.status = "completed";
    if any(string({cycle_receipts.status}) == "failed")
        receipt.status = "completed_with_errors";
    end
    receipt.completed_sections = completed_section_names(execution_plan);
    receipt.cycle_receipts = cycle_receipts;
    receipt.reference_roi_file = reference_roi_file;
    receipt.completed_cycle_count = sum(string({cycle_receipts.status}) == "completed");
    receipt.failed_cycle_count = sum(string({cycle_receipts.status}) == "failed");
    receipt.skipped_cycle_count = sum(string({cycle_receipts.status}) == "skipped");

    record_average_receipt = struct();
    if request.record.record_average
        average_control = struct( ...
            'mode', "ap", ...
            'record_path', string(request.input_path), ...
            'cycle_results', cycle_receipts, ...
            'reference_roi_file', reference_roi_file, ...
            'source_results_path', string(request.workflow.source_results_path), ...
            'output_dir_name', "AP_analysis3_rec_average", ...
            'request', request, ...
            'run_frequency', request.params.ap.record_frequency.enabled, ...
            'frequency', rmfield(request.params.ap.record_frequency,'enabled'), ...
            'discovery_max_depth', request.workflow.discovery_max_depth, ...
            'discovery_max_mat_files', request.workflow.discovery_max_mat_files);
        record_average_receipt = Record_analysis3(average_control, runtime);
        receipt.completed_sections(end+1,1) = "record_average";
    end
    receipt.record_average = record_average_receipt;

    summary_file = fullfile(summary_path, 'AP_analysis3_rec_batch_summary.mat');
    output_files = string(summary_file);
    if strlength(diary_file) > 0
        output_files(end + 1, 1) = diary_file;
    end
    if isstruct(record_average_receipt) && isfield(record_average_receipt,'output_files')
        output_files = unique([output_files;string(record_average_receipt.output_files(:))],'stable');
    end
    receipt.output_files = output_files;
    record_request = request; %#ok<NASGU>
    record_execution_plan = execution_plan; %#ok<NASGU>
    save(summary_file, 'receipt', 'cycle_receipts', 'record_average_receipt', ...
        'record_request', 'record_execution_plan', 'items', '-v7.3');
    update_manifest_if_present(manifest_file, struct( ...
        'status', receipt.status, ...
        'completed_sections', receipt.completed_sections, ...
        'failed_section', "", ...
        'output_files', output_files, ...
        'completed_at', datetime('now'), ...
        'error', struct(), ...
        'cycle_receipts', cycle_receipts, ...
        'reference_roi_file', reference_roi_file, ...
        'record_average', record_average_receipt));
    emit_event_safe(runtime, make_event( ...
        "record_completed", request, summary_path, ...
        sprintf('AP Record scheduling completed (%d completed, %d failed, %d skipped).', ...
        receipt.completed_cycle_count, receipt.failed_cycle_count, receipt.skipped_cycle_count)));
catch ME
    receipt.status = "failed";
    receipt.error = exception_to_struct(ME);
    receipt.cycle_receipts = cycle_receipts;
    receipt.failed_cycle_count = sum(string({cycle_receipts.status}) == "failed");
    update_manifest_if_present(manifest_file, struct( ...
        'status', "failed", ...
        'failed_section', "cycle_dispatch", ...
        'completed_at', datetime('now'), ...
        'error', receipt.error, ...
        'cycle_receipts', cycle_receipts));
    emit_event_safe(runtime, make_event( ...
        "record_failed", request, summary_path, ME.message));
    rethrow(ME);
end
end

function cycle_receipt = run_one_cycle(cycle_request, runtime, item)
emit_event_safe(runtime, struct( ...
    'timestamp', datetime('now'), ...
    'type', "cycle_started", ...
    'mode', "ap", ...
    'scope', "record", ...
    'section', "cycle_dispatch", ...
    'state', "running", ...
    'message', "Starting " + string(item.label) + ...
        "; reuse_source=" + string(cycle_request.workflow.source_results_path), ...
    'input_path', string(item.cycle_path), ...
    'output_path', string(cycle_request.workflow.output_path)));
cycle_receipt = aaa.workflows.AP_analysis3(cycle_request, runtime);
cycle_receipt.label = string(item.label);
cycle_receipt.cycle_name = string(item.cycle_name);
cycle_receipt.record_path = string(item.record_path);
end

function cycle_request = build_cycle_request(record_request, item, run_tag, roi_file,runtime)
cycle_request = record_request;
cycle_request.scope = "single";
cycle_request.input_path = string(item.cycle_path);
cycle_request.record = struct();
cycle_run_name = string(item.record_name) + "_" + ...
    string(item.cycle_name) + "_" + string(run_tag);
cycle_request.workflow.run_name = cycle_run_name;

configured_root = string(record_request.workflow.output_path);
if strlength(configured_root) > 0
    cycle_output_root = fullfile(configured_root, char(item.label));
else
    cycle_output_root = fullfile(item.cycle_path, char(record_request.record.output_root));
end
cycle_request.workflow.output_path = string(fullfile(cycle_output_root, char(cycle_run_name)));

if strlength(roi_file) > 0
    cycle_request.workflow.sections.roi = "reuse";
    cycle_request.workflow.roi_path = string(roi_file);
end
source_root = string(record_request.workflow.source_results_path);
cycle_request.workflow.source_results_path = source_root;
[~,match_kind,preflight_cycle]=aaa.helpers.common.runtime_preflight(runtime,cycle_request);
if match_kind=="cycle"&&strlength(string(preflight_cycle.selected_reuse_source))>0
    cycle_request.workflow.source_results_path=string(preflight_cycle.selected_reuse_source);
end
[cycle_request, resolution] = aaa.helpers.common.resolve_request_reuse_source( ...
    cycle_request, struct('cycle_name',item.cycle_name,'label',item.label));
if resolution.attempted && strlength(resolution.source_path) == 0
    error('AAA:APRecord:ReuseSourceMissing', ...
        'No complete AP result satisfies the reuse plan for %s under %s.', ...
        item.label, resolution.report.search_root);
end
end

function roi_file = resolve_initial_reference_roi(request, items, execution_plan)
roi_file = "";
if execution_plan.roi ~= "run" || request.record.redraw_reference_roi
    return;
end
explicit_path = string(request.record.reference_roi_path);
if strlength(explicit_path) == 0
    explicit_path = string(request.workflow.roi_path);
end
if strlength(explicit_path) > 0
    roi_file = find_roi_file(explicit_path);
    return;
end
if ~request.record.reuse_reference_roi
    return;
end
reference_idx = find(string({items.cycle_name}) == ...
    string(request.record.reference_cycle_name), 1, 'first');
if isempty(reference_idx)
    return;
end
candidates = [dir(fullfile(items(reference_idx).cycle_path,'**','results','roi.mat')); ...
    dir(fullfile(items(reference_idx).cycle_path,'**','1_raw_ROI.mat'))];
if ~isempty(candidates)
    [~, newest_idx] = max([candidates.datenum]);
    roi_file = string(fullfile(candidates(newest_idx).folder, candidates(newest_idx).name));
end
end

function roi_file = find_roi_file(path_value)
path_value = string(path_value);
if isfile(path_value)
    roi_file = path_value;
    return;
end
if ~isfolder(path_value)
    error('AAA:APRecord:ROINotFound', 'Reference ROI path does not exist: %s', path_value);
end
preferred=[string(fullfile(path_value,'results','roi.mat')); ...
    string(fullfile(path_value,'1_raw_ROI.mat'))];
for idx=1:numel(preferred)
    if isfile(preferred(idx)),roi_file=preferred(idx);return;end
end
candidates = [dir(fullfile(path_value,'**','results','roi.mat')); ...
    dir(fullfile(path_value,'**','1_raw_ROI.mat'))];
if isempty(candidates)
    error('AAA:APRecord:ROINotFound', ...
        'No unified or legacy ROI result found under %s.', path_value);
end
[~, newest_idx] = max([candidates.datenum]);
roi_file = string(fullfile(candidates(newest_idx).folder, candidates(newest_idx).name));
end

function receipt = maybe_existing_receipt(cycle_request, item)
receipt = empty_cycle_receipt();
explicit_file = fullfile(cycle_request.workflow.output_path, '-1_explicit_results.mat');
manifest_file = fullfile(cycle_request.workflow.output_path, 'analysis_manifest.mat');
has_result_marker = isfile(explicit_file) || isfile(manifest_file);
if has_result_marker && aaa.io.manifest_allows_reuse(cycle_request.workflow.output_path)
    receipt.status = "skipped";
    receipt.mode = "ap";
    receipt.scope = "single";
    receipt.input_path = string(cycle_request.input_path);
    receipt.output_path = string(cycle_request.workflow.output_path);
    receipt.manifest_file = string(manifest_file);
    receipt.label = string(item.label);
    receipt.cycle_name = string(item.cycle_name);
    receipt.record_path = string(item.record_path);
end
end

function receipt = failed_cycle_receipt(item, cycle_request, ME)
receipt = empty_cycle_receipt();
receipt.status = "failed";
receipt.mode = "ap";
receipt.scope = "single";
receipt.input_path = string(cycle_request.input_path);
receipt.output_path = string(cycle_request.workflow.output_path);
receipt.label = string(item.label);
receipt.cycle_name = string(item.cycle_name);
receipt.record_path = string(item.record_path);
receipt.failed_section = "legacy_ap_workflow";
receipt.error = exception_to_struct(ME);
end

function receipt = empty_record_receipt(request, output_path, manifest_file)
receipt = struct( ...
    'status', "running", ...
    'mode', "ap", ...
    'scope', "record", ...
    'input_path', string(request.input_path), ...
    'output_path', string(output_path), ...
    'manifest_file', string(manifest_file), ...
    'completed_sections', strings(0, 1), ...
    'failed_section', "", ...
    'error', struct(), ...
    'cycle_receipts', repmat(empty_cycle_receipt(), 0, 1), ...
    'reference_roi_file', "", ...
    'completed_cycle_count', 0, ...
    'failed_cycle_count', 0, ...
    'skipped_cycle_count', 0, ...
    'record_average', struct(), ...
    'output_files', strings(0, 1));
end

function receipt = empty_cycle_receipt()
receipt = struct( ...
    'receipt_schema_version', "1.0.0", ...
    'status', "pending", ...
    'mode', "ap", ...
    'scope', "single", ...
    'input_path', "", ...
    'output_path', "", ...
    'started_at', NaT, ...
    'completed_at', NaT, ...
    'completed_sections', strings(0, 1), ...
    'reused_sections', strings(0, 1), ...
    'skipped_sections', strings(0, 1), ...
    'failed_section', "", ...
    'failed_channel', "", ...
    'error', struct(), ...
    'section_status', struct(), ...
    'channel_status', struct([]), ...
    'manifest_file', "", ...
    'roi_file', "", ...
    'output_files', strings(0, 1), ...
    'label', "", ...
    'cycle_name', "", ...
    'record_path', "");
end

function event = make_event(event_type, request, output_path, message)
state = canonical_event_state(event_type);
event = struct( ...
    'timestamp', datetime('now'), ...
    'type', string(event_type), ...
    'mode', "ap", ...
    'scope', "record", ...
    'section', "record", ...
    'state', state, ...
    'message', string(message), ...
    'input_path', string(request.input_path), ...
    'output_path', string(output_path));
end

function state = canonical_event_state(event_type)
event_type = lower(string(event_type));
if contains(event_type, "failed")
    state = "failed";
elseif contains(event_type, "completed")
    state = "completed";
elseif contains(event_type, "started")
    state = "running";
else
    state = "";
end
end

function emit_event_safe(runtime, event)
try
    aaa.io.emit_event(runtime, event);
catch ME
    warning('AAA:APRecord:EventDeliveryFailed', ...
        'AP Record event delivery failed: %s', ME.message);
end
end

function failure = exception_to_struct(ME)
failure = struct( ...
    'identifier', string(ME.identifier), ...
    'message', string(ME.message), ...
    'report', string(getReport(ME, 'extended', 'hyperlinks', 'off')));
end

function names = completed_section_names(plan)
all_names = string(fieldnames(plan));
keep = false(size(all_names));
for idx = 1:numel(all_names)
    keep(idx) = string(plan.(char(all_names(idx)))) ~= "skip";
end
names = all_names(keep);
end

function close_owned_logger(logger, owns_logger)
if owns_logger && ~isempty(logger) && isvalid(logger)
    logger.close();
end
end

function update_manifest_if_present(manifest_file, patch)
if strlength(string(manifest_file)) > 0 && isfile(manifest_file)
    aaa.io.update_manifest(manifest_file, patch);
end
end

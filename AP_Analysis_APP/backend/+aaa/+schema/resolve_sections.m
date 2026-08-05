function [canonical_preset, preset_sections, overlay_sections, resolved_sections, resolution_info] = ...
    resolve_sections(mode, requested_preset, overlay_sections)
%RESOLVE_SECTIONS Resolve an AP or Dual preset into one validated plan.

mode = normalize_mode(mode);
catalog = aaa.schema.section_catalog(mode);
section_names = cellstr([catalog.name]);
valid_actions = ["run", "reuse", "skip"];

if nargin < 2 || isempty(requested_preset)
    requested_preset = "full";
end
if nargin < 3 || isempty(overlay_sections)
    overlay_sections = struct();
end
if ~isstruct(overlay_sections) || ~isscalar(overlay_sections)
    error('AAA:Schema:InvalidSectionOverlay', ...
        'Section overlay must be one scalar struct.');
end

requested_preset = lower(strtrim(string(requested_preset)));
if ~isscalar(requested_preset) || strlength(requested_preset) == 0
    error('AAA:Schema:InvalidPreset', ...
        'Preset must be one nonempty scalar string.');
end

[canonical_preset, legacy_alias] = canonicalize_preset(requested_preset);
preset_sections = preset_plan(mode, canonical_preset);

overlay_names = fieldnames(overlay_sections);
unknown_names = setdiff(overlay_names, section_names);
if ~isempty(unknown_names)
    error('AAA:Schema:UnknownSection', ...
        'Unknown %s section field(s): %s', mode, strjoin(unknown_names, ', '));
end
overlay_sections = normalize_actions(overlay_sections, overlay_names, valid_actions);

if canonical_preset == "custom"
    missing_names = setdiff(section_names, overlay_names);
    if ~isempty(missing_names)
        error('AAA:Schema:IncompleteCustomSections', ...
            ['preset="custom" requires every section action explicitly. ' ...
             'Missing field(s): %s'], strjoin(missing_names, ', '));
    end
    resolved_sections = reorder_sections(overlay_sections, section_names);
else
    resolved_sections = preset_sections;
    for idx = 1:numel(overlay_names)
        name = overlay_names{idx};
        resolved_sections.(name) = overlay_sections.(name);
    end
end

validate_dependencies(mode, resolved_sections);
resolution_info = struct( ...
    'mode', mode, ...
    'requested_preset', requested_preset, ...
    'canonical_preset', canonical_preset, ...
    'legacy_alias', legacy_alias, ...
    'custom_requires_complete_overlay', canonical_preset == "custom", ...
    'section_names', string(section_names));
end

function [canonical, legacy_alias] = canonicalize_preset(requested)
canonical = requested;
legacy_alias = "";
switch requested
    case "custom1"
        canonical = "reuse_motion";
        legacy_alias = requested;
    case "custom2"
        canonical = "reuse_trace";
        legacy_alias = requested;
    case "custom3"
        canonical = "skip_motion";
        legacy_alias = requested;
end
if strlength(legacy_alias) > 0
    warning('AAA:Schema:LegacyPresetAlias', ...
        'Preset "%s" is deprecated; use "%s".', legacy_alias, canonical);
end
end

function plan = preset_plan(mode, preset)
switch mode
    case "dual"
        switch preset
            case "full"
                plan = dual_plan("run","run","run","run","run","run","run","run","run","run","run");
            case "reuse_motion"
                plan = dual_plan("run","reuse","skip","run","run","run","run","run","run","skip","run");
            case "reuse_trace"
                plan = dual_plan("skip","skip","skip","skip","reuse","run","skip","skip","run","skip","run");
            case "skip_motion"
                plan = dual_plan("run","skip","skip","run","run","run","skip","skip","run","skip","run");
            case "analysis_only"
                plan = dual_plan("skip","skip","skip","skip","reuse","reuse","run","run","run","run","run");
            case "motion_only"
                plan = dual_plan("run","run","skip","skip","skip","skip","skip","skip","skip","skip","skip");
            case "custom"
                plan = struct();
            otherwise
                error('AAA:Schema:UnknownPreset', ...
                    ['Unknown Dual preset: %s. Use full, reuse_motion, reuse_trace, ' ...
                     'skip_motion, analysis_only, motion_only, or custom.'], preset);
        end
    case "ap"
        switch preset
            case "full"
                plan = ap_plan("run","run","run","run","run","run","skip","skip","run","run");
            case "reuse_motion"
                plan = ap_plan("run","reuse","run","run","run","run","skip","skip","run","run");
            case "reuse_trace"
                plan = ap_plan("skip","skip","skip","reuse","run","run","skip","skip","run","run");
            case "skip_motion"
                plan = ap_plan("run","skip","run","run","run","run","skip","skip","run","run");
            case "analysis_only"
                plan = ap_plan("skip","skip","skip","reuse","reuse","run","skip","skip","run","run");
            case "motion_only"
                plan = ap_plan("run","run","skip","skip","skip","skip","skip","skip","skip","skip");
            case "custom"
                plan = struct();
            otherwise
                error('AAA:Schema:UnknownPreset', ...
                    ['Unknown AP preset: %s. Use full, reuse_motion, reuse_trace, ' ...
                     'skip_motion, analysis_only, motion_only, or custom.'], preset);
        end
end
end

function plan = dual_plan(input,motion,registration,roi,trace,peak,visualization,comparison,stim,frequency,export)
plan = struct('input',input,'motion',motion,'registration',registration, ...
    'roi',roi,'trace',trace,'peak',peak,'visualization',visualization, ...
    'comparison',comparison,'stim',stim,'frequency',frequency,'export',export);
end

function plan = ap_plan(input,motion,roi,trace,peak,stim,ap_events,ap_statistics,visualization,export)
plan = struct('input',input,'motion',motion,'roi',roi,'trace',trace, ...
    'peak',peak,'stim',stim,'ap_events',ap_events,'ap_statistics',ap_statistics, ...
    'visualization',visualization,'export',export);
end

function sections = normalize_actions(sections, names, valid_actions)
for idx = 1:numel(names)
    name = names{idx};
    action = lower(strtrim(string(sections.(name))));
    if ~isscalar(action) || ~ismember(action, valid_actions)
        error('AAA:Schema:InvalidSectionAction', ...
            'Section %s has invalid action %s. Use run, reuse, or skip.', name, action);
    end
    sections.(name) = action;
end
end

function sections = reorder_sections(input_sections, names)
sections = struct();
for idx = 1:numel(names)
    name = names{idx};
    sections.(name) = input_sections.(name);
end
end

function validate_dependencies(mode, sections)
if sections.export == "reuse"
    error('AAA:Schema:InvalidSectionAction', ...
        'export="reuse" is unsupported. Use run or skip.');
end
if sections.input == "reuse"
    error('AAA:Schema:InvalidSectionAction', ...
        'input="reuse" is unsupported; saved result restoration belongs to trace="reuse".');
end
if sections.input == "skip"
    upstream=[string(sections.motion),string(sections.roi)];
    if mode=="dual",upstream(end+1)=string(sections.registration);end
    if any(upstream ~= "skip") || ~ismember(string(sections.trace),["reuse","skip"])
        error('AAA:Schema:InvalidSectionDependency', ...
            'input="skip" requires motion/registration/roi="skip" and trace="reuse" or "skip".');
    end
end
if ismember(sections.motion, ["run", "reuse"]) && sections.input ~= "run"
    error('AAA:Schema:InvalidSectionDependency', ...
        'motion="run" or "reuse" requires input="run".');
end
if sections.roi == "run" && sections.input ~= "run"
    error('AAA:Schema:InvalidSectionDependency', ...
        'roi="run" requires input="run".');
end
roi_available = ismember(sections.roi, ["run", "reuse"]);
if sections.trace == "run" && ~(sections.input == "run" && roi_available)
    error('AAA:Schema:InvalidSectionDependency', ...
        'trace="run" requires input="run" and roi="run" or "reuse".');
end
trace_available = ismember(sections.trace, ["run", "reuse"]);
if sections.peak ~= "skip" && ~trace_available
    error('AAA:Schema:InvalidSectionDependency', ...
        'peak requires trace="run" or "reuse".');
end

if mode == "dual"
    if sections.registration == "run" && sections.input ~= "run"
        error('AAA:Schema:InvalidSectionDependency', ...
            'registration="run" requires input="run".');
    end
    if sections.registration == "reuse" && sections.input ~= "run"
        error('AAA:Schema:InvalidSectionDependency', ...
            'registration="reuse" requires input="run".');
    end
    invalidate_reused_downstream(sections, 'motion', {'trace'});
    invalidate_reused_downstream(sections, 'registration', {'roi','trace'});
    invalidate_reused_downstream(sections, 'roi', {'trace'});
    invalidate_reused_downstream(sections, 'trace', {'peak','visualization','comparison','stim','frequency'});
    downstream = ["visualization","comparison","stim","frequency"];
else
    invalidate_reused_downstream(sections, 'motion', {'trace'});
    invalidate_reused_downstream(sections, 'roi', {'trace'});
    invalidate_reused_downstream(sections, 'trace', {'peak','stim','ap_events','ap_statistics','visualization'});
    invalidate_reused_downstream(sections, 'peak', {'ap_events','ap_statistics'});
    downstream = ["stim","ap_events","ap_statistics","visualization"];
    peak_available = ismember(sections.peak, ["run","reuse"]);
    if sections.ap_events ~= "skip" && ~peak_available
        error('AAA:Schema:InvalidSectionDependency', ...
            'ap_events requires peak="run" or "reuse".');
    end
    events_available = ismember(sections.ap_events, ["run","reuse"]);
    if sections.ap_statistics ~= "skip" && ~events_available
        error('AAA:Schema:InvalidSectionDependency', ...
            'ap_statistics requires ap_events="run" or "reuse".');
    end
end
for name = downstream
    if sections.(name) == "run" && ~trace_available
        error('AAA:Schema:InvalidSectionDependency', ...
            '%s="run" requires trace="run" or "reuse".', name);
    end
end

end

function invalidate_reused_downstream(sections, upstream, downstream_names)
if sections.(upstream) ~= "run"
    return;
end
for idx = 1:numel(downstream_names)
    downstream = downstream_names{idx};
    if sections.(downstream) == "reuse"
        error('AAA:Schema:InvalidSectionDependency', ...
            ['Invalid plan: %s="run" invalidates %s="reuse". ' ...
             'Set %s to run or skip.'], upstream, downstream, downstream);
    end
end
end

function mode = normalize_mode(mode)
mode = lower(strtrim(string(mode)));
if ~isscalar(mode) || ~ismember(mode, ["ap", "dual"])
    error('AAA:Schema:InvalidMode', 'mode must be "ap" or "dual".');
end
end

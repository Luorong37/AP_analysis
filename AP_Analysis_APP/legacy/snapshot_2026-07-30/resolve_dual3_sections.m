function [canonical_preset, preset_sections, overlay_sections, resolved_sections, resolution_info] = ...
    resolve_dual3_sections(requested_preset, overlay_sections)
%RESOLVE_DUAL3_SECTIONS Resolve one Dual_analysis3 preset into an execution plan.
%
% Configuration layers:
%   preset_sections   complete named-preset baseline (empty for custom)
%   overlay_sections  caller-supplied section overrides, preserved as given
%   resolved_sections complete validated execution plan used by the pipeline

section_names = { ...
    'input', 'motion', 'registration', 'roi', 'trace', 'peak', ...
    'visualization', 'comparison', 'stim', 'frequency', 'export'};
valid_actions = ["run", "reuse", "skip"];

if nargin < 1 || isempty(requested_preset)
    requested_preset = "full";
end
if nargin < 2 || isempty(overlay_sections)
    overlay_sections = struct();
end
if ~isstruct(overlay_sections) || ~isscalar(overlay_sections)
    error('Dual_analysis3:InvalidSectionOverlay', ...
        'Section overlay must be one scalar struct.');
end

requested_preset = lower(strtrim(string(requested_preset)));
if ~isscalar(requested_preset) || strlength(requested_preset) == 0
    error('Dual_analysis3:InvalidPreset', ...
        'Preset must resolve to one nonempty scalar string.');
end

canonical_preset = requested_preset;
legacy_alias = "";
switch requested_preset
    case "custom1"
        canonical_preset = "reuse_motion";
        legacy_alias = requested_preset;
    case "custom2"
        canonical_preset = "reuse_trace";
        legacy_alias = requested_preset;
    case "custom3"
        canonical_preset = "skip_motion";
        legacy_alias = requested_preset;
end
if strlength(legacy_alias) > 0
    warning('Dual_analysis3:LegacyPresetAlias', ...
        'Preset "%s" is deprecated; use "%s".', legacy_alias, canonical_preset);
end

switch canonical_preset
    case "full"
        preset_sections = make_sections( ...
            "run", "run", "run", "run", "run", "run", ...
            "run", "run", "run", "run", "run");
    case "reuse_motion"
        preset_sections = make_sections( ...
            "run", "reuse", "skip", "run", "run", "run", ...
            "run", "run", "run", "skip", "run");
    case "reuse_trace"
        preset_sections = make_sections( ...
            "reuse", "skip", "skip", "reuse", "reuse", "run", ...
            "skip", "skip", "run", "skip", "run");
    case "skip_motion"
        preset_sections = make_sections( ...
            "run", "skip", "skip", "run", "run", "run", ...
            "skip", "skip", "run", "skip", "run");
    case "analysis_only"
        preset_sections = make_sections( ...
            "reuse", "skip", "skip", "skip", "reuse", "reuse", ...
            "run", "run", "run", "run", "run");
    case "motion_only"
        preset_sections = make_sections( ...
            "run", "run", "skip", "skip", "skip", "skip", ...
            "skip", "skip", "skip", "skip", "skip");
    case "custom"
        preset_sections = struct();
    otherwise
        error('Dual_analysis3:UnknownPreset', ...
            ['Unknown preset: %s. Use full, reuse_motion, reuse_trace, ' ...
             'skip_motion, analysis_only, motion_only, or custom.'], canonical_preset);
end

overlay_names = fieldnames(overlay_sections);
unknown_names = setdiff(overlay_names, section_names);
if ~isempty(unknown_names)
    error('Dual_analysis3:UnknownSection', ...
        'Unknown section overlay field(s): %s', strjoin(unknown_names, ', '));
end
overlay_sections = normalize_section_actions(overlay_sections, overlay_names, valid_actions);

if canonical_preset == "custom"
    missing_names = setdiff(section_names, overlay_names);
    if ~isempty(missing_names)
        error('Dual_analysis3:IncompleteCustomSections', ...
            ['preset="custom" requires all section actions explicitly. ' ...
             'Missing field(s): %s'], strjoin(missing_names, ', '));
    end
    resolved_sections = reorder_sections(overlay_sections, section_names);
else
    resolved_sections = preset_sections;
    for idx = 1:numel(overlay_names)
        section_name = overlay_names{idx};
        resolved_sections.(section_name) = overlay_sections.(section_name);
    end
end

validate_section_dependencies(resolved_sections);
resolution_info = struct( ...
    'requested_preset', requested_preset, ...
    'canonical_preset', canonical_preset, ...
    'legacy_alias', legacy_alias, ...
    'custom_requires_complete_overlay', canonical_preset == "custom", ...
    'section_names', {string(section_names)});
end

function sections = make_sections(input, motion, registration, roi, trace, peak, ...
    visualization, comparison, stim, frequency, export)
sections = struct( ...
    'input', input, ...
    'motion', motion, ...
    'registration', registration, ...
    'roi', roi, ...
    'trace', trace, ...
    'peak', peak, ...
    'visualization', visualization, ...
    'comparison', comparison, ...
    'stim', stim, ...
    'frequency', frequency, ...
    'export', export);
end

function sections = normalize_section_actions(sections, section_names, valid_actions)
for idx = 1:numel(section_names)
    section_name = section_names{idx};
    action = lower(strtrim(string(sections.(section_name))));
    if ~isscalar(action) || ~ismember(action, valid_actions)
        error('Dual_analysis3:InvalidSectionAction', ...
            'Section %s has invalid action %s. Use run, reuse, or skip.', ...
            section_name, action);
    end
    sections.(section_name) = action;
end
end

function sections = reorder_sections(sections_in, section_names)
sections = struct();
for idx = 1:numel(section_names)
    section_name = section_names{idx};
    sections.(section_name) = sections_in.(section_name);
end
end

function validate_section_dependencies(sections)
dependency_edges = {
    'motion',       {'trace'};
    'registration', {'roi', 'trace'};
    'roi',          {'trace'};
    'trace',        {'peak', 'visualization', 'comparison', 'stim', 'frequency'}};
for edge_idx = 1:size(dependency_edges, 1)
    upstream = dependency_edges{edge_idx, 1};
    downstream_list = dependency_edges{edge_idx, 2};
    if sections.(upstream) == "run"
        for downstream_idx = 1:numel(downstream_list)
            downstream = downstream_list{downstream_idx};
            if sections.(downstream) == "reuse"
                error('Dual_analysis3:InvalidSectionDependency', ...
                    ['Invalid execution plan: %s="run" invalidates %s="reuse". ' ...
                     'Set %s to "run" or "skip".'], upstream, downstream, downstream);
            end
        end
    end
end

trace_available = ismember(sections.trace, ["run", "reuse"]);
roi_available = ismember(sections.roi, ["run", "reuse"]);
if ismember(sections.motion, ["run", "reuse"]) && sections.input ~= "run"
    error('Dual_analysis3:InvalidSectionDependency', ...
        'motion="run" or "reuse" requires input="run".');
end
if sections.registration == "run" && sections.input ~= "run"
    error('Dual_analysis3:InvalidSectionDependency', ...
        'registration="run" requires input="run".');
end
if sections.registration == "reuse" && sections.input ~= "reuse"
    error('Dual_analysis3:InvalidSectionDependency', ...
        ['registration="reuse" requires input="reuse". ' ...
         'Use registration="skip" when reloading raw inputs.']);
end
if sections.roi == "run" && sections.input ~= "run"
    error('Dual_analysis3:InvalidSectionDependency', ...
        'roi="run" requires input="run".');
end
if sections.trace == "run" && ~(sections.input == "run" && roi_available)
    error('Dual_analysis3:InvalidSectionDependency', ...
        'trace="run" requires input="run" and roi="run" or "reuse".');
end
for source_section = ["trace", "peak", "visualization", "comparison", "stim", "frequency"]
    if sections.(source_section) == "reuse" && sections.input ~= "reuse"
        error('Dual_analysis3:InvalidSectionDependency', ...
            '%s="reuse" requires input="reuse" so source containers are initialized atomically.', ...
            source_section);
    end
end
if sections.peak ~= "skip" && ~trace_available
    error('Dual_analysis3:InvalidSectionDependency', ...
        'peak requires trace="run" or "reuse".');
end
for dependent_name = ["visualization", "comparison", "stim", "frequency"]
    if sections.(dependent_name) == "run" && ~trace_available
        error('Dual_analysis3:InvalidSectionDependency', ...
            '%s="run" requires trace="run" or "reuse".', dependent_name);
    end
end
if sections.input == "skip"
    non_input_actions = rmfield(sections, 'input');
    if any(structfun(@(value) string(value) ~= "skip", non_input_actions))
        error('Dual_analysis3:InvalidSectionDependency', ...
            'input="skip" is only legal when every other section is also "skip".');
    end
end
if sections.export == "reuse"
    error('Dual_analysis3:InvalidSectionAction', ...
        'export="reuse" is unsupported. Use export="run" or "skip".');
end
end

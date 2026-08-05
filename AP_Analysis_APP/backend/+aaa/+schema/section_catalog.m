function catalog = section_catalog(mode)
%SECTION_CATALOG Return the ordered AAA workflow-stage definitions.
% The catalog is the single source used by the App, resolver, validator,
% workflow status table, and manifest section-status fields.

if nargin < 1 || isempty(mode)
    mode = "dual";
end
mode = normalize_mode(mode);

common_actions = ["run", "reuse", "skip"];
export_actions = ["run", "skip"];

switch mode
    case "dual"
        rows = {
            'input',         'Input',                  'Load raw/rebuilt movies; saved trace containers belong to Trace reuse.';
            'motion',        'Motion correction',      'Estimate, apply, reuse, or skip motion shifts.';
            'registration',  'Channel registration',   'Register the calcium and voltage image coordinates.';
            'roi',           'ROI',                    'Create/select or reuse ROI definitions.';
            'trace',         'Trace processing',       'Extract and process voltage/calcium traces.';
            'peak',          'Voltage peaks',          'Detect, edit, or reuse accepted voltage peaks.';
            'stim',          'Stimulus analysis',      'Build stimulus windows and stimulus-aware results.';
            'visualization', 'Visualization',          'Generate non-essential summary figures and images.';
            'comparison',    'Dual comparison',        'Generate voltage-calcium comparisons.';
            'frequency',     'Frequency analysis',     'Generate FFT and wavelet results.';
            'export',        'Export',                 'Write the explicit final result bundle.'};
    case "ap"
        rows = {
            'input',          'Input',              'Load a single voltage movie; saved AP traces belong to Trace reuse.';
            'motion',         'Motion correction',  'Estimate, apply, reuse, or skip motion shifts.';
            'roi',            'ROI',                'Create/select or reuse ROI definitions.';
            'trace',          'Trace processing',   'Extract and process AP voltage traces.';
            'peak',           'Peak detection',     'Detect, edit, or reuse accepted AP peaks.';
            'stim',           'Stimulus analysis',  'Build voltage stimulus windows and stimulus-aware results.';
            'ap_events',      'AP events (pending)', 'Excluded in this migration phase; reserved for event windows/FWHM.';
            'ap_statistics',  'AP statistics (pending)', 'Excluded in this migration phase; reserved for summary and ISI analyses.';
            'visualization',  'Visualization',      'Generate non-essential summary figures and images.';
            'export',         'Export',             'Write the explicit final result bundle.'};
end

template = struct( ...
    'name', "", ...
    'mode', mode, ...
    'order', 0, ...
    'ui_label', "", ...
    'description', "", ...
    'allowed_actions', strings(0, 1), ...
    'default_action_control', "dropdown", ...
    'is_visualization_only', false);
catalog = repmat(template, size(rows, 1), 1);
for idx = 1:size(rows, 1)
    catalog(idx) = template;
    catalog(idx).name = string(rows{idx, 1});
    catalog(idx).order = idx;
    catalog(idx).ui_label = string(rows{idx, 2});
    catalog(idx).description = string(rows{idx, 3});
    if mode == "ap" && ismember(catalog(idx).name,["ap_events","ap_statistics"])
        catalog(idx).allowed_actions = "skip";
        catalog(idx).default_action_control = "readonly";
    elseif catalog(idx).name == "export"
        catalog(idx).allowed_actions = export_actions;
    elseif catalog(idx).name == "input"
        catalog(idx).allowed_actions = ["run","skip"];
    else
        catalog(idx).allowed_actions = common_actions;
    end
    catalog(idx).is_visualization_only = catalog(idx).name == "visualization";
end
end

function mode = normalize_mode(mode)
mode = lower(strtrim(string(mode)));
if ~isscalar(mode) || ~ismember(mode, ["ap", "dual"])
    error('AAA:Schema:InvalidMode', 'mode must be "ap" or "dual".');
end
end

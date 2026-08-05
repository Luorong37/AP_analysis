function [ctx, outcome] = roi(ctx, action)
%ROI Coordinate shared AP/Dual ROI selection, reuse, and persistence.
% Scientific selection and trace extraction live in ordinary select_rois.m.

action = normalize_action(action);
roles = channel_roles(ctx);
outcome = struct('section',"roi",'action',action,'state',"completed", ...
    'channels',roles(:),'message',"");
if action == "skip"
    outcome.state = "skipped";
    outcome.message = "ROI section skipped.";
    return;
end
validate_context(ctx);
output_path = require_output_path(ctx);

% Load inputs: obtain channel movies/profiles or reusable ROI masks.
[movies,profiles] = channel_inputs(ctx);

if action == "reuse"
    aaa.sections.emit_phase(ctx,"roi","load_masks","running", ...
        "Loading reusable ROI masks.");
    [saved_spec,source_file,compatibility_report] = load_roi_spec(ctx);
    offset_xy = resolve_offset(ctx,saved_spec);
    spec = struct('masks',saved_spec.masks,'offset_xy',offset_xy);
    map_results = struct();
    selection_source = "saved ROI masks";
else
    offset_xy = resolve_offset(ctx,struct());
    aaa.sections.emit_phase(ctx,"roi","activity_maps","running", ...
        "Creating role-aware activity maps.");
    [maps,map_results] = create_channel_maps(movies,profiles,roles);
    spec = struct('maps',maps,'offset_xy',offset_xy);
    spec = apply_runtime_selection(spec,ctx.runtime);
    source_file = "";
    selection_source = "new ROI selection";
    compatibility_report = struct();
end

% Run algorithm: apply supplied masks or perform the configured selection.
interactive=action~="reuse"&&(~isstruct(ctx.runtime)|| ...
    ~isfield(ctx.runtime,'roi_spec')||isempty(ctx.runtime.roi_spec));
if interactive
    aaa.sections.emit_phase(ctx,"roi","manual_selection","waiting", ...
        "Waiting for manual ROI selection.");
else
    aaa.sections.emit_phase(ctx,"roi","apply_masks","running", ...
        "Applying supplied or reusable ROI masks.");
end
[rois,traces,selection_info] = select_rois(movies,profiles,spec);
aaa.sections.emit_phase(ctx,"roi","extract_traces","completed", ...
    "ROI masks and raw traces are available.");
selection_info.source_file = string(source_file);
selection_info.selection_source = selection_source;
nrois = validate_trace_pairing(traces,roles);

% Return results: register masks, raw traces, maps, and provenance.
roi_file = string(fullfile(output_path,'results','roi.mat'));
for idx = 1:numel(ctx.channels)
    role = char(roles(idx));
    mask = mask_for_role(rois,roles(idx));
    ctx.channels(idx).data.roi_mask = mask;
    ctx.channels(idx).results = store_trace_stage( ...
        ctx.channels(idx).results,'raw',traces.(role),{},roi_file, ...
        ctx.channels(idx).results.movie_info,'select_rois',struct( ...
            'role',roles(idx),'offset_xy',offset_xy, ...
            'selection_mode',selection_info.selection_mode));
    ctx.channels(idx).output_files.roi = strings(0,1);
end
ctx.shared.results.roi = struct('rois',rois,'info',selection_info, ...
    'offset_xy',offset_xy,'nrois',nrois);
if ~isempty(fieldnames(map_results))
    ctx.shared.results.maps = map_results;
end

aaa.sections.emit_phase(ctx,"roi","persist","running", ...
    "Preparing the unified ROI result.");
outcome.compatibility_report=compatibility_report;
outcome.message = sprintf('%s completed with %d paired ROI(s).', ...
    selection_source,nrois);
end

function action = normalize_action(action)
action = lower(strtrim(string(action)));
if ~isscalar(action) || ~ismember(action,["run","reuse","skip"])
    error('AAA:Sections:InvalidRoiAction', ...
        'ROI action must be run, reuse, or skip.');
end
end

function validate_context(ctx)
if ~isstruct(ctx) || ~isscalar(ctx) || ~isfield(ctx,'channels') ...
        || isempty(ctx.channels)
    error('AAA:Sections:InvalidRoiContext', ...
        'ROI context must contain one or two channels.');
end
for idx = 1:numel(ctx.channels)
    channel = ctx.channels(idx);
    if ~isfield(channel,'profile') || ~isfield(channel.profile,'role') ...
            || ~isfield(channel.profile,'roi') ...
            || ~isfield(channel.data,'movie_3d') ...
            || isempty(channel.data.movie_3d) ...
            || ~isfield(channel.results,'movie_info')
        error('AAA:Sections:IncompleteRoiChannel', ...
            'ROI requires loaded movie, profile, and movie_info for channel %d.',idx);
    end
end
end

function roles = channel_roles(ctx)
roles = strings(0,1);
if ~isstruct(ctx) || ~isfield(ctx,'channels') || isempty(ctx.channels), return; end
profiles = [ctx.channels.profile];
roles = string({profiles.role});
end

function [movies,profiles] = channel_inputs(ctx)
movies = cell(1,numel(ctx.channels));
for idx = 1:numel(ctx.channels), movies{idx} = ctx.channels(idx).data.movie_3d; end
profiles = [ctx.channels.profile];
end

function [maps,results] = create_channel_maps(movies,profiles,roles)
maps = struct(); results = struct();
for idx = 1:numel(roles)
    role = char(roles(idx));
    [maps.(role),map_info] = create_map(movies{idx},profiles(idx));
    results.(role) = struct('data',maps.(role),'info',map_info);
end
end

function spec = apply_runtime_selection(spec,runtime)
% Runtime injection is restricted to selection geometry for App/tests.
% Role, profiles, maps, and offset remain section-owned single authorities.
if ~isstruct(runtime) || ~isfield(runtime,'roi_spec') ...
        || isempty(runtime.roi_spec)
    return;
end
value = runtime.roi_spec;
if ~isstruct(value) || ~isscalar(value)
    error('AAA:Sections:InvalidRuntimeRoiSpec', ...
        'runtime.roi_spec must be one scalar struct.');
end
allowed = {'polygons','masks'};
unknown = setdiff(fieldnames(value),allowed);
if ~isempty(unknown)
    error('AAA:Sections:InvalidRuntimeRoiSpec', ...
        ['runtime.roi_spec may supply only polygons or masks. Role, maps, ' ...
         'and offset are bound by the request/context.']);
end
for idx = 1:numel(allowed)
    if isfield(value,allowed{idx}), spec.(allowed{idx}) = value.(allowed{idx}); end
end
end

function offset = resolve_offset(ctx,saved_spec)
offset = [];
if isfield(ctx,'shared') && isfield(ctx.shared,'results') ...
        && isfield(ctx.shared.results,'registration') ...
        && isfield(ctx.shared.results.registration,'offset_xy')
    offset = ctx.shared.results.registration.offset_xy;
end
if isempty(offset) && isfield(ctx.request.workflow,'reuse_offset')
    offset = ctx.request.workflow.reuse_offset;
end
if isempty(offset) && isfield(saved_spec,'offset_xy')
    offset = saved_spec.offset_xy;
end
if isempty(offset), offset = [0 0]; end
if ~isnumeric(offset) || numel(offset) ~= 2 || any(~isfinite(offset))
    error('AAA:Sections:InvalidRoiOffset', ...
        'Resolved ROI offset must be one finite [x y] vector.');
end
offset = double(reshape(offset,1,2));
end

function [spec,source_file,report] = load_roi_spec(ctx)
source = string(ctx.request.workflow.roi_path);
if strlength(source) == 0, source = string(ctx.request.workflow.source_results_path); end
if strlength(source) == 0
    error('AAA:Sections:RoiReuseSourceMissing', ...
        'roi="reuse" requires workflow.roi_path or source_results_path.');
end
[results,report]=load_analysis_results(source,struct( ...
    'mode',ctx.mode,'scope',ctx.scope,'runtime',ctx.runtime));
if ~isfield(results,'roi') || ~isstruct(results.roi)
    error('AAA:Sections:RoiReuseMissing', ...
        'No compatible ROI result was found at %s.',source);
end
saved=results.roi;
if isfield(saved,'rois'),saved_for_masks=struct('rois',saved.rois);else,saved_for_masks=saved;end
[masks,found]=extract_saved_masks(saved_for_masks,string(ctx.mode));
if ~found
    error('AAA:Sections:RoiReuseMissing', ...
        'The saved ROI result has no compatible masks: %s.',source);
end
spec=struct('masks',masks,'offset_xy',[]);
if isfield(saved,'offset_xy'),spec.offset_xy=saved.offset_xy;
elseif isfield(saved,'offset'),spec.offset_xy=saved.offset;end
source_file="";
if ~isempty(report.absolute_source_files)
    source_file=report.absolute_source_files(1);
end
end

function [masks,found] = extract_saved_masks(saved,mode)
masks = struct(); found = false;
if isfield(saved,'rois') && isstruct(saved.rois)
    if isfield(saved.rois,'bwmask') && ~isempty(saved.rois.bwmask)
        masks.voltage = saved.rois.bwmask; found = true;
    end
    if isfield(saved.rois,'bwmask_ca') && ~isempty(saved.rois.bwmask_ca)
        masks.calcium = saved.rois.bwmask_ca; found = true;
    end
end
if ~found && isfield(saved,'bwmask') && ~isempty(saved.bwmask)
    masks.voltage = saved.bwmask; found = true;
elseif ~found && isfield(saved,'mask') && ~isempty(saved.mask)
    masks.voltage = saved.mask; found = true;
end
if mode == "ap" && isfield(masks,'calcium') && ~isfield(masks,'voltage')
    masks.voltage = masks.calcium;
    masks = rmfield(masks,'calcium');
end
end

function nrois = validate_trace_pairing(traces,roles)
nrois = [];
for role = roles
    count = size(traces.(char(role)),2);
    if isempty(nrois), nrois = count;
    elseif count ~= nrois
        error('AAA:Sections:RoiTracePairingMismatch', ...
            'All channels must produce the same number of ROI trace columns.');
    end
end
if isempty(nrois) || nrois < 1
    error('AAA:Sections:EmptyRoiResult','ROI section produced no ROI traces.');
end
end

function mask = mask_for_role(rois,role)
if string(role) == "calcium"
    mask = rois.bwmask_ca;
else
    mask = rois.bwmask;
end
end

function output_path = require_output_path(ctx)
output_path = string(ctx.output_path);
if strlength(output_path) == 0
    error('AAA:Sections:MissingOutputPath','ROI requires a nonempty output path.');
end
if ~isfolder(output_path), mkdir(output_path); end
end

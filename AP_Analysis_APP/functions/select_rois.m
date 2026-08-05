function [rois, traces, info] = select_rois(movies, profiles, spec)
%SELECT_ROIS Select or replay labeled ROIs for one or two channels.
% This is the shared AP/Dual ROI functional unit. Channel behavior is read
% from profile.role. In dual mode the frozen coordinate convention is:
%   voltage_position = calcium_position + offset_xy.
%
% MOVIES is a cell array paired with the profile struct array. SPEC fields:
%   maps       struct keyed by role (optional guidance images)
%   masks      struct keyed by role (replay route)
%   polygons   struct array with role and N-by-2 position [x y]
%   offset_xy  [x y], default [0 0]
% If neither masks nor polygons are supplied, MATLAB Desktop opens the
% interactive selector. Returned TRACES is keyed by role and preserves each
% movie's own frame count; no temporal alignment is performed.
% Input formats:
%   movies   1-by-C cell array (C=1 AP, C=2 Dual); each movie is accepted
%            by normalize_movie_matrix and all channels share frame_size.
%   profiles 1-by-C bound profile struct array with unique role fields.
%   spec.masks fields are [ncols x nrows] nonnegative integer label maps;
%            polygon positions are N-by-2 [x y] coordinates.
% Outputs:
%   rois    Struct keyed by role containing 2-D integer label maps.
%   traces  Struct keyed by role containing [frames x ROI] matrices.
%   info    Scalar selection/provenance struct.

if nargin < 3 || isempty(spec), spec = struct(); end
[movies, profiles, roles, frame_size] = normalize_inputs(movies, profiles);
spec = normalize_spec(spec, roles);

masks = initialize_masks(spec, roles, frame_size);
if has_any_mask(masks, roles)
    masks = reconcile_masks(masks, roles, spec.offset_xy, frame_size);
    selection_mode = "reused_masks";
elseif ~isempty(spec.polygons)
    masks = masks_from_polygons(spec.polygons, roles, ...
        spec.offset_xy, frame_size);
    selection_mode = "supplied_polygons";
else
    masks = interactive_masks(movies, profiles, roles, ...
        spec.maps, spec.offset_xy, frame_size);
    selection_mode = "interactive_polygons";
end

[rois, traces, extraction] = build_outputs( ...
    movies, profiles, roles, masks, frame_size);
info = struct( ...
    'roles', roles, ...
    'selection_mode', selection_mode, ...
    'offset_xy', spec.offset_xy, ...
    'coordinate_rule', "voltage_position = calcium_position + offset_xy", ...
    'frame_size', frame_size, ...
    'trace_alignment_rule', "none; each channel keeps its own frame count", ...
    'extraction', extraction, ...
    'created_at', datetime("now"));
end

function [movies, profiles, roles, frame_size] = normalize_inputs(movies, profiles)
if ~iscell(movies) || isempty(movies)
    error('AAA:Functions:InvalidRoiMovies', ...
        'movies must be a nonempty cell array.');
end
if ~isstruct(profiles) || numel(profiles) ~= numel(movies) ...
        || ~ismember(numel(movies), [1 2])
    error('AAA:Functions:InvalidRoiProfiles', ...
        'profiles must contain one entry per movie (one AP or two Dual).');
end
profiles = reshape(profiles,1,[]);
movies = reshape(movies,1,[]);
roles = strings(1,numel(profiles));
frame_size = [];
for idx = 1:numel(profiles)
    [~, movie_info] = normalize_movie_matrix(movies{idx},profiles(idx));
    roles(idx) = lower(string(profiles(idx).role));
    if isempty(frame_size)
        frame_size = movie_info.frame_size;
    elseif ~isequal(frame_size,movie_info.frame_size)
        error('AAA:Functions:RoiGeometryMismatch', ...
            'All ROI channels must have the same analysis frame size.');
    end
end
if numel(unique(roles)) ~= numel(roles) ...
        || (numel(roles) == 2 && ~isequal(sort(roles),["calcium","voltage"]))
    error('AAA:Functions:InvalidRoiRoles', ...
        'Dual ROI requires one voltage and one calcium profile.');
end
end

function spec = normalize_spec(spec, roles)
if ~isstruct(spec) || ~isscalar(spec)
    error('AAA:Functions:InvalidRoiSpec', ...
        'ROI spec must be one scalar struct.');
end
defaults = struct('maps',struct(),'masks',struct(), ...
    'polygons',repmat(struct('role',"",'position',[]),0,1), ...
    'offset_xy',[0 0]);
fields = fieldnames(defaults);
for idx = 1:numel(fields)
    if ~isfield(spec,fields{idx}), spec.(fields{idx}) = defaults.(fields{idx}); end
end
if ~isstruct(spec.maps) || ~isscalar(spec.maps) ...
        || ~isstruct(spec.masks) || ~isscalar(spec.masks)
    error('AAA:Functions:InvalidRoiSpec', ...
        'spec.maps and spec.masks must be scalar structs keyed by role.');
end
if ~isnumeric(spec.offset_xy) || numel(spec.offset_xy) ~= 2 ...
        || any(~isfinite(spec.offset_xy))
    error('AAA:Functions:InvalidRoiOffset', ...
        'spec.offset_xy must be a finite [x y] vector.');
end
spec.offset_xy = double(reshape(spec.offset_xy,1,2));
if ~isstruct(spec.polygons)
    error('AAA:Functions:InvalidRoiPolygons', ...
        'spec.polygons must be a struct array.');
end
for idx = 1:numel(spec.polygons)
    polygon = spec.polygons(idx);
    if ~all(isfield(polygon,{'role','position'})) ...
            || ~ismember(lower(string(polygon.role)),roles) ...
            || ~isnumeric(polygon.position) ...
            || size(polygon.position,2) ~= 2 ...
            || size(polygon.position,1) < 3 ...
            || any(~isfinite(polygon.position),'all')
        error('AAA:Functions:InvalidRoiPolygons', ...
            'Every polygon needs a valid role and at least three finite [x y] vertices.');
    end
end
end

function masks = initialize_masks(spec, roles, frame_size)
masks = struct();
for role = roles
    name = char(role);
    if isfield(spec.masks,name) && ~isempty(spec.masks.(name))
        mask = double(spec.masks.(name));
        validate_mask(mask,frame_size,role);
        masks.(name) = mask;
    else
        masks.(name) = zeros(frame_size);
    end
end
end

function tf = has_any_mask(masks, roles)
tf = false;
for role = roles
    tf = tf || any(masks.(char(role)) > 0,'all');
end
end

function masks = reconcile_masks(masks, roles, offset_xy, frame_size)
if isscalar(roles)
    require_labels(masks.(char(roles)),roles);
    return;
end
has_voltage = any(masks.voltage > 0,'all');
has_calcium = any(masks.calcium > 0,'all');
if has_voltage && ~has_calcium
    masks.calcium = translate_label_mask( ...
        masks.voltage,-offset_xy(1),-offset_xy(2));
elseif has_calcium && ~has_voltage
    masks.voltage = translate_label_mask( ...
        masks.calcium,offset_xy(1),offset_xy(2));
end
validate_mask(masks.voltage,frame_size,"voltage");
validate_mask(masks.calcium,frame_size,"calcium");
labels_voltage = positive_labels(masks.voltage);
labels_calcium = positive_labels(masks.calcium);
if isempty(labels_voltage) || ~isequal(labels_voltage,labels_calcium)
    error('AAA:Functions:RoiPairingMismatch', ...
        ['Voltage and calcium masks must contain the same nonempty labels. ' ...
         'A translated ROI may have fallen completely outside one image.']);
end
end

function masks = masks_from_polygons(polygons, roles, offset_xy, frame_size)
masks = struct();
for role = roles, masks.(char(role)) = zeros(frame_size); end
for roi_idx = 1:numel(polygons)
    source_role = lower(string(polygons(roi_idx).role));
    positions = positions_for_roles( ...
        polygons(roi_idx).position,source_role,roles,offset_xy);
    for role = roles
        name = char(role);
        mask = polygon_mask(positions.(name),frame_size);
        if ~any(mask,'all')
            error('AAA:Functions:RoiOutsideImage', ...
                'ROI %d does not overlap the %s image.',roi_idx,role);
        end
        masks.(name)(mask) = roi_idx;
    end
end
end

function positions = positions_for_roles(position, source_role, roles, offset_xy)
positions = struct();
positions.(char(source_role)) = double(position);
if numel(roles) == 1, return; end
if source_role == "voltage"
    positions.calcium = translate_position(position,-offset_xy);
else
    positions.voltage = translate_position(position,offset_xy);
end
end

function position = translate_position(position, offset_xy)
position = double(position);
position(:,1) = position(:,1) + offset_xy(1);
position(:,2) = position(:,2) + offset_xy(2);
end

function mask = polygon_mask(position, frame_size)
position = double(position);
position(:,1) = max(1,min(frame_size(2),position(:,1)));
position(:,2) = max(1,min(frame_size(1),position(:,2)));
mask = poly2mask(position(:,1),position(:,2),frame_size(1),frame_size(2));
end

function masks = interactive_masks(movies, profiles, roles, maps, offset_xy, frame_size)
if ~usejava('desktop')
    error('AAA:Functions:InteractiveRoiRequiresDesktop', ...
        ['Interactive ROI selection requires MATLAB Desktop. Supply ' ...
         'spec.masks or spec.polygons for a noninteractive call.']);
end
means = struct();
matrices = cell(1,numel(roles));
for idx = 1:numel(roles)
    [matrices{idx},~] = normalize_movie_matrix(movies{idx},profiles(idx));
    means.(char(roles(idx))) = reshape( ...
        mean(double(matrices{idx}),2),frame_size);
end
masks = struct();
for role = roles, masks.(char(role)) = zeros(frame_size); end

fig = figure('Color','w','Name','AAA ROI Selection', ...
    'Position',get(0,'Screensize'));
cleanup = onCleanup(@() close_if_valid(fig)); %#ok<NASGU>
[image_axes, trace_axes] = build_gui(fig,roles,maps,means);
colors = lines(100);
while isgraphics(fig)
    sgtitle(fig,['Click a map/image and draw a polygon. ' ...
        'Enter or Q: finish; R: remove last ROI.']);
    set(fig,'CurrentCharacter',char(0));
    waitforbuttonpress;
    key = lower(get(fig,'CurrentCharacter'));
    if any(key == [char(13),'q']), break; end
    if key == 'r'
        masks = remove_last_masks(masks,roles);
        refresh_overlays(image_axes,masks,roles,colors);
        continue;
    end
    [source_role,valid_axis] = selected_role(gca,image_axes,roles);
    if ~valid_axis, continue; end
    roi_idx = max(all_labels(masks,roles),[],'all') + 1;
    polygon = drawpolygon('Parent',gca,'Color', ...
        colors(mod(roi_idx-1,size(colors,1))+1,:));
    position = double(polygon.Position);
    delete(polygon);
    if size(position,1) < 3, continue; end
    positions = positions_for_roles(position,source_role,roles,offset_xy);
    for role = roles
        mask = polygon_mask(positions.(char(role)),frame_size);
        masks.(char(role))(mask) = roi_idx;
    end
    refresh_overlays(image_axes,masks,roles,colors);
    refresh_traces(trace_axes,matrices,masks,roles,roi_idx,colors);
end
if isempty(positive_labels(masks.(char(roles(1)))))
    error('AAA:Functions:NoRoisSelected','No ROI was selected.');
end
masks = reconcile_masks(masks,roles,offset_xy,frame_size);
end

function [image_axes,trace_axes] = build_gui(fig,roles,maps,means)
image_axes = struct();
trace_axes = struct();
layout = tiledlayout(fig,numel(roles),3,'TileSpacing','compact');
for role = roles
    name = char(role);
    map_ax = nexttile(layout);
    if isfield(maps,name) && ~isempty(maps.(name))
        imagesc(map_ax,maps.(name)); axis(map_ax,'image'); axis(map_ax,'ij');
    else
        text(map_ax,.5,.5,'Map unavailable','Units','normalized', ...
            'HorizontalAlignment','center'); axis(map_ax,'off');
    end
    title(map_ax,role + " map"); hold(map_ax,'on');
    image_ax = nexttile(layout);
    imagesc(image_ax,means.(name)); axis(image_ax,'image'); axis(image_ax,'ij');
    title(image_ax,role + " mean"); hold(image_ax,'on');
    trace_ax = nexttile(layout); grid(trace_ax,'on'); title(trace_ax,role + " trace");
    image_axes.(name) = [map_ax,image_ax];
    trace_axes.(name) = trace_ax;
end
end

function [role,valid] = selected_role(ax,image_axes,roles)
role = ""; valid = false;
for candidate = roles
    if any(ax == image_axes.(char(candidate)))
        role = candidate; valid = true; return;
    end
end
end

function masks = remove_last_masks(masks,roles)
labels = all_labels(masks,roles);
if isempty(labels), return; end
last = max(labels);
for role = roles
    name = char(role);
    masks.(name)(masks.(name) == last) = 0;
end
end

function refresh_overlays(image_axes,masks,roles,colors)
for role = roles
    name = char(role);
    for ax = image_axes.(name)
        delete(findall(ax,'Tag','AAA_ROI_OVERLAY'));
        for label = positive_labels(masks.(name))
            boundaries = bwboundaries(masks.(name) == label);
            if isempty(boundaries), continue; end
            boundary = boundaries{1};
            line(ax,boundary(:,2),boundary(:,1),'Tag','AAA_ROI_OVERLAY', ...
                'Color',colors(mod(label-1,size(colors,1))+1,:));
        end
    end
end
end

function refresh_traces(trace_axes,matrices,masks,roles,roi_idx,colors)
for idx = 1:numel(roles)
    role = roles(idx); name = char(role);
    selected = masks.(name)(:) == roi_idx;
    trace = mean(matrices{idx}(selected,:),1,'double');
    cla(trace_axes.(name));
    plot(trace_axes.(name),trace,'Color', ...
        colors(mod(roi_idx-1,size(colors,1))+1,:));
    title(trace_axes.(name),sprintf('%s ROI %d',role,roi_idx));
end
end

function [rois,traces,extraction] = build_outputs( ...
        movies,profiles,roles,masks,frame_size)
traces = struct(); extraction = struct();
for idx = 1:numel(roles)
    role = char(roles(idx));
    [traces.(role),extraction.(role)] = extract_roi_traces( ...
        movies{idx},masks.(role),profiles(idx));
end
labels = positive_labels(masks.(char(roles(1))));
for role = roles
    if ~isequal(labels,positive_labels(masks.(char(role))))
        error('AAA:Functions:RoiPairingMismatch', ...
            'All channel masks must contain identical labels.');
    end
end
rois = struct('bwmask',masks.(char(roles(1))), ...
    'boundary',{{}},'position',{{}});
if ismember("voltage",roles), rois.bwmask = masks.voltage; end
if ismember("calcium",roles), rois.bwmask_ca = masks.calcium; end
for label = labels
    boundary = first_boundary(rois.bwmask == label);
    rois.boundary{end+1} = boundary; %#ok<AGROW>
    rois.position{end+1} = [boundary(:,2),boundary(:,1)]; %#ok<AGROW>
end
if ismember("calcium",roles)
    rois.boundary_ca = {};
    rois.position_ca = {};
    for label = labels
        boundary = first_boundary(rois.bwmask_ca == label);
        rois.boundary_ca{end+1} = boundary; %#ok<AGROW>
        rois.position_ca{end+1} = [boundary(:,2),boundary(:,1)]; %#ok<AGROW>
    end
end
if ~isequal(size(rois.bwmask),frame_size)
    error('AAA:Functions:RoiOutputGeometryMismatch', ...
        'Internal ROI output geometry mismatch.');
end
end

function boundary = first_boundary(mask)
boundaries = bwboundaries(mask);
if isempty(boundaries)
    error('AAA:Functions:EmptyRoiLabel','An ROI label contains no pixels.');
end
boundary = boundaries{1};
end

function translated = translate_label_mask(mask,dx,dy)
if exist('imtranslate','file') == 2
    translated = imtranslate(double(mask),[dx dy],'nearest','FillValues',0);
else
    translated = zeros(size(mask));
    sx = round(dx); sy = round(dy);
    rows = max(1,1-sy):min(size(mask,1),size(mask,1)-sy);
    cols = max(1,1-sx):min(size(mask,2),size(mask,2)-sx);
    if ~isempty(rows) && ~isempty(cols)
        translated(rows+sy,cols+sx) = mask(rows,cols);
    end
end
translated = round(translated);
translated(translated < 0) = 0;
end

function validate_mask(mask,frame_size,role)
if ~isnumeric(mask) || ~isequal(size(mask),frame_size) ...
        || any(~isfinite(mask),'all') || any(mask < 0,'all') ...
        || any(mask ~= round(mask),'all')
    error('AAA:Functions:InvalidRoiMask', ...
        '%s mask must be a [%d %d] nonnegative integer label matrix.', ...
        role,frame_size(1),frame_size(2));
end
end

function require_labels(mask,role)
if isempty(positive_labels(mask))
    error('AAA:Functions:EmptyRoiMask','%s mask has no positive labels.',role);
end
end

function labels = positive_labels(mask)
labels = unique(double(mask(:)),'sorted')';
labels = labels(labels > 0);
end

function labels = all_labels(masks,roles)
labels = [];
for role = roles
    labels = union(labels,positive_labels(masks.(char(role))));
end
end

function close_if_valid(fig)
if isgraphics(fig), close(fig); end
end

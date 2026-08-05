function [rois, traces, traces_ca, offset] = select_ROI_dual(movie, movie_ca, ...
    nrows, ncols, correct, map, map_ca, mask, mask_ca, offset)

% select_ROI_dual
% Interactive paired ROI selection for dual-channel voltage/calcium data.
%
% Inputs keep the historical Dual_analysis3 interface:
%   movie, movie_ca : [pixels x frames]
%   nrows, ncols    : analysis geometry passed in from Dual_analysis3
%   correct         : offset estimation mode
%                     'none'            -> keep reuse_offset / [0 0]
%                     'manual_points'   -> click one matching point on each
%                                          channel average image
%                     'matlab_register' -> use MATLAB built-in translation
%                                          registration on average images
%                     true/false are still accepted for legacy callers
%   map, map_ca     : optional guidance maps
%   mask, mask_ca   : optional pre-existing paired ROI masks
%   offset          : optional [x_offset, y_offset], where
%                     voltage_position = calcium_position + offset
%
% Default behavior is intentionally conservative:
%   if no valid offset is provided and correct == 'none', the function uses
%   zero offset and does not force the user into manual correction.
%
% Offset convention:
%   voltage_position = calcium_position + offset
% So any estimated offset is the shift that moves calcium into the voltage
% coordinate system.

if nargin < 5 || isempty(correct)
    correct = false;
end
if nargin < 6
    map = [];
end
if nargin < 7
    map_ca = [];
end
if nargin < 8
    mask = [];
end
if nargin < 9
    mask_ca = [];
end
if nargin < 10
    offset = [];
end

colors = lines(100);

[image_height, image_width] = resolve_dual_image_size(movie, nrows, ncols);
% Keep full movies in their input type; only the small averages and ROI
% traces need double output.
mean_voltage = reshape(mean(movie, 2, 'double'), image_height, image_width);
mean_calcium = reshape(mean(movie_ca, 2, 'double'), image_height, image_width);

offset_mode = normalize_dual_offset_mode(correct);
offset = sanitize_dual_offset(offset);
x_offset = offset(1);
y_offset = offset(2);

traces = [];
traces_ca = [];
rois = struct( ...
    'bwmask', zeros(image_height, image_width), ...
    'bwmask_ca', zeros(image_height, image_width), ...
    'boundary', {{}}, ...
    'boundary_ca', {{}}, ...
    'position', {{}}, ...
    'position_ca', {{}});

fig = figure('Color', 'w', 'Name', 'Dual ROI Selection');
set(fig, 'Position', get(0, 'Screensize'));
set(fig, 'KeyPressFcn', @(src, event) set_dual_key_pressed(event, fig));

axes_info = create_dual_roi_gui(fig, map, map_ca, mean_voltage, mean_calcium);
update_dual_offset_instruction(axes_info, offset_mode, x_offset, y_offset, ...
    'Click voltage or calcium map/image to draw ROI');

% Only estimate a new offset automatically when explicitly requested by the
% selected offset mode.
if ~strcmpi(offset_mode, "none")
    [x_offset, y_offset] = estimate_dual_offset( ...
        axes_info, mean_voltage, mean_calcium, image_height, image_width, offset_mode);
    offset = [x_offset, y_offset];
    update_dual_offset_instruction(axes_info, offset_mode, x_offset, y_offset, ...
        'Offset estimated. Click voltage or calcium map/image to draw ROI');
else
    offset = [x_offset, y_offset];
end

% Reuse an existing ROI set after the offset decision is finalized.
% This ordering matters: when the user asked for manual_points or
% matlab_register, we should estimate the new offset first and then rebuild
% the calcium-side mask from the voltage ROI geometry instead of silently
% replaying the old paired masks.
if ~isempty(mask) || ~isempty(mask_ca)
    [mask, mask_ca, reuse_note] = reconcile_reused_dual_masks(mask, mask_ca, x_offset, y_offset, offset_mode);
    if strlength(reuse_note) > 0
        disp(char(reuse_note));
    end
    if ~isempty(mask) && ~isempty(mask_ca)
        [rois, traces, traces_ca] = replay_dual_masks( ...
            mask, mask_ca, movie, movie_ca, axes_info, colors, image_height, image_width);
        offset = [x_offset, y_offset];
        disp('Finished ROIs selection');
        return;
    end
end

overlay_handles = cell(0, 1);

while true
    set(fig, 'CurrentCharacter', char(0));
    waitforbuttonpress;

    pressed_char = get(fig, 'CurrentCharacter');
    if ~isempty(pressed_char)
        switch lower(pressed_char)
            case {char(13), 'q'}
                disp('Finished ROIs selection');
                return;
            case 'r'
                [rois, traces, traces_ca, overlay_handles] = remove_last_dual_roi( ...
                    rois, traces, traces_ca, overlay_handles, axes_info);
                continue;
            case 'o'
                if strcmpi(offset_mode, "none")
                    update_dual_offset_instruction(axes_info, offset_mode, x_offset, y_offset, ...
                        'Offset mode is none. Change correct_offset_mode in Dual_analysis3 to estimate a new offset.');
                    continue;
                end
                [x_offset, y_offset] = estimate_dual_offset( ...
                    axes_info, mean_voltage, mean_calcium, image_height, image_width, offset_mode);
                offset = [x_offset, y_offset];
                update_dual_offset_instruction(axes_info, offset_mode, x_offset, y_offset, ...
                    'Offset updated. Click voltage or calcium map/image to draw ROI');
                continue;
        end
    end

    current_axis = gca;
    if ~is_dual_roi_axis(current_axis, axes_info)
        continue;
    end

    num_roi = max(rois.bwmask(:)) + 1;
    color = colors(mod(num_roi - 1, size(colors, 1)) + 1, :);

    [selected_channel, mask_v, mask_ca, boundary_v, boundary_ca, position_v, position_ca] = ...
        collect_dual_roi_pair(current_axis, axes_info, image_height, image_width, x_offset, y_offset, color);

    rois.bwmask(mask_v) = num_roi;
    rois.bwmask_ca(mask_ca) = num_roi;
    rois.boundary{end+1} = boundary_v; %#ok<AGROW>
    rois.boundary_ca{end+1} = boundary_ca; %#ok<AGROW>
    rois.position{end+1} = position_v; %#ok<AGROW>
    rois.position_ca{end+1} = position_ca; %#ok<AGROW>

    trace_v = mean(movie(mask_v(:), :), 1, 'double')';
    trace_ca = mean(movie_ca(mask_ca(:), :), 1, 'double')';
    traces(:, end+1) = trace_v; %#ok<AGROW>
    traces_ca(:, end+1) = trace_ca; %#ok<AGROW>

    overlay_handles{end+1, 1} = plot_dual_roi_overlay( ... %#ok<AGROW>
        axes_info, boundary_v, boundary_ca, num_roi, color);

    cla(axes_info.trace_v);
    plot(axes_info.trace_v, trace_v, 'Color', color, 'LineWidth', 1);
    xlabel(axes_info.trace_v, 'Frames');
    ylabel(axes_info.trace_v, 'Voltage Intensity');
    grid(axes_info.trace_v, 'on');
    title(axes_info.trace_v, sprintf('Voltage ROI %d (%s-selected)', num_roi, selected_channel));

    cla(axes_info.trace_ca);
    plot(axes_info.trace_ca, trace_ca, 'Color', color, 'LineWidth', 1);
    xlabel(axes_info.trace_ca, 'Frames');
    ylabel(axes_info.trace_ca, 'Calcium Intensity');
    grid(axes_info.trace_ca, 'on');
    title(axes_info.trace_ca, sprintf('Calcium ROI %d', num_roi));

    update_dual_offset_instruction(axes_info, offset_mode, x_offset, y_offset, ...
        sprintf('ROI %d saved. Continue drawing or finish.', num_roi));

    key = wait_for_dual_key(fig);
    switch lower(string(key))
        case {"return", "q"}
            disp('Finished ROIs selection');
            return;
        case {"space", ""}
            continue;
        case "r"
            [rois, traces, traces_ca, overlay_handles] = remove_last_dual_roi( ...
                rois, traces, traces_ca, overlay_handles, axes_info);
            continue;
        case "o"
            if strcmpi(offset_mode, "none")
                update_dual_offset_instruction(axes_info, offset_mode, x_offset, y_offset, ...
                    'Offset mode is none. Change correct_offset_mode in Dual_analysis3 to estimate a new offset.');
                continue;
            end
            [x_offset, y_offset] = estimate_dual_offset( ...
                axes_info, mean_voltage, mean_calcium, image_height, image_width, offset_mode);
            offset = [x_offset, y_offset];
            update_dual_offset_instruction(axes_info, offset_mode, x_offset, y_offset, ...
                'Offset updated. Continue drawing or finish.');
            continue;
        otherwise
            continue;
    end
end

end

function [image_height, image_width] = resolve_dual_image_size(movie, nrows, ncols)
expected_pixels = nrows * ncols;
if size(movie, 1) ~= expected_pixels
    error('Movie size (%d pixels) does not match nrows*ncols (%d).', size(movie, 1), expected_pixels);
end

% Dual_analysis3 historically stores movies as [ncols x nrows x frames].
% Keep that geometry convention here so ROI masks map back to the movie
% exactly the same way as the rest of the pipeline.
image_height = ncols;
image_width = nrows;
end

function offset = sanitize_dual_offset(offset)
if isempty(offset)
    offset = [0, 0];
elseif numel(offset) ~= 2 || any(~isfinite(offset))
    error('Offset must be empty or a finite 1x2 vector [x_offset, y_offset].');
else
    offset = reshape(double(offset), 1, 2);
end
end

function axes_info = create_dual_roi_gui(fig, map, map_ca, mean_voltage, mean_calcium)
clf(fig);

axes_info = struct();

axes_info.map_v = subplot(2, 3, 1, 'Parent', fig);
display_dual_panel(axes_info.map_v, map, 'Voltage Map');

axes_info.image_v = subplot(2, 3, 2, 'Parent', fig);
display_dual_panel(axes_info.image_v, mean_voltage, 'Voltage Image');

axes_info.trace_v = subplot(2, 3, 3, 'Parent', fig);
hold(axes_info.trace_v, 'on');
grid(axes_info.trace_v, 'on');
xlabel(axes_info.trace_v, 'Frames');
ylabel(axes_info.trace_v, 'Voltage Intensity');
title(axes_info.trace_v, 'Voltage Trace');

axes_info.map_ca = subplot(2, 3, 4, 'Parent', fig);
display_dual_panel(axes_info.map_ca, map_ca, 'Calcium Map');

axes_info.image_ca = subplot(2, 3, 5, 'Parent', fig);
display_dual_panel(axes_info.image_ca, mean_calcium, 'Calcium Image');

axes_info.trace_ca = subplot(2, 3, 6, 'Parent', fig);
hold(axes_info.trace_ca, 'on');
grid(axes_info.trace_ca, 'on');
xlabel(axes_info.trace_ca, 'Frames');
ylabel(axes_info.trace_ca, 'Calcium Intensity');
title(axes_info.trace_ca, 'Calcium Trace');
end

function display_dual_panel(ax, image_data, panel_title)
axes(ax); %#ok<LAXES>
cla(ax);
if isempty(image_data)
    axis(ax, 'off');
    text(ax, 0.5, 0.5, 'Unavailable', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    title(ax, panel_title);
    return;
end

imagesc(ax, image_data);
axis(ax, 'image');
axis(ax, 'ij');
colormap(ax, gray);
title(ax, panel_title);
hold(ax, 'on');
end

function update_dual_instruction(ax, message_text)
title(ax, message_text, 'Interpreter', 'none');
end

function update_dual_offset_instruction(axes_info, offset_mode, x_offset, y_offset, status_text)
if strcmpi(offset_mode, "none")
    offset_action = 'O: disabled in none mode';
else
    offset_action = 'O: re-estimate offset';
end

update_dual_instruction(axes_info.trace_v, sprintf([ ...
    'Dual ROI selection\n' ...
    '%s\n' ...
    'SPACE: next ROI | R: remove last ROI | ENTER/Q: finish | %s\n' ...
    'Offset mode = %s | Current offset [x y] = [%.2f %.2f]'], ...
    status_text, offset_action, char(offset_mode), x_offset, y_offset));
end

function tf = is_dual_roi_axis(ax, axes_info)
tf = isequal(ax, axes_info.map_v) ...
    || isequal(ax, axes_info.map_ca) ...
    || isequal(ax, axes_info.image_v) ...
    || isequal(ax, axes_info.image_ca);
end

function [selected_channel, mask_v, mask_ca, boundary_v, boundary_ca, position_v, position_ca] = ...
    collect_dual_roi_pair(current_axis, axes_info, image_height, image_width, x_offset, y_offset, color)

is_voltage_axis = isequal(current_axis, axes_info.map_v) || isequal(current_axis, axes_info.image_v);
if is_voltage_axis
    selected_channel = 'voltage';
    position_v = draw_dual_polygon(current_axis, color);
    position_ca = translate_dual_polygon(position_v, -x_offset, -y_offset);
else
    selected_channel = 'calcium';
    position_ca = draw_dual_polygon(current_axis, color);
    position_v = translate_dual_polygon(position_ca, x_offset, y_offset);
end

[mask_v, boundary_v, position_v] = polygon_to_dual_mask(position_v, image_height, image_width);
[mask_ca, boundary_ca, position_ca] = polygon_to_dual_mask(position_ca, image_height, image_width);
end

function position = draw_dual_polygon(ax, color)
while true
    roi_polygon = drawpolygon('Color', color, 'LineWidth', 1.5, 'Parent', ax);
    position = double(roi_polygon.Position);
    delete(roi_polygon);
    if size(position, 1) >= 3
        return;
    end
    uiwait(msgbox('Polygon must have at least 3 vertices. Please draw again.', ...
        'Invalid ROI', 'warn'));
end
end

function position_out = translate_dual_polygon(position_in, dx, dy)
position_out = double(position_in);
position_out(:, 1) = position_out(:, 1) + dx;
position_out(:, 2) = position_out(:, 2) + dy;
end

function [mask, boundary, position] = polygon_to_dual_mask(position, image_height, image_width)
position = double(position);
position(:, 1) = max(1, min(image_width, position(:, 1)));
position(:, 2) = max(1, min(image_height, position(:, 2)));

mask = poly2mask(position(:, 1), position(:, 2), image_height, image_width);
boundary_cells = bwboundaries(mask);
if isempty(boundary_cells)
    error('The ROI polygon does not overlap the image after offset mapping.');
end
boundary = boundary_cells{1};
end

function handles = plot_dual_roi_overlay(axes_info, boundary_v, boundary_ca, roi_index, color)
handles = gobjects(8, 1);

handles(1) = plot(axes_info.map_v, boundary_v(:, 2), boundary_v(:, 1), 'Color', color, 'LineWidth', 1.2);
handles(2) = plot(axes_info.image_v, boundary_v(:, 2), boundary_v(:, 1), 'Color', color, 'LineWidth', 1.2);
handles(3) = plot(axes_info.map_ca, boundary_ca(:, 2), boundary_ca(:, 1), 'Color', color, 'LineWidth', 1.2);
handles(4) = plot(axes_info.image_ca, boundary_ca(:, 2), boundary_ca(:, 1), 'Color', color, 'LineWidth', 1.2);

handles(5) = text(axes_info.map_v, mean(boundary_v(:, 2)), mean(boundary_v(:, 1)), num2str(roi_index), ...
    'Color', color, 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
handles(6) = text(axes_info.image_v, mean(boundary_v(:, 2)), mean(boundary_v(:, 1)), num2str(roi_index), ...
    'Color', color, 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
handles(7) = text(axes_info.map_ca, mean(boundary_ca(:, 2)), mean(boundary_ca(:, 1)), num2str(roi_index), ...
    'Color', color, 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
handles(8) = text(axes_info.image_ca, mean(boundary_ca(:, 2)), mean(boundary_ca(:, 1)), num2str(roi_index), ...
    'Color', color, 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

function [rois, traces, traces_ca] = replay_dual_masks(mask, mask_ca, movie, movie_ca, axes_info, colors, image_height, image_width)
mask = double(mask);
mask_ca = double(mask_ca);

rois = struct( ...
    'bwmask', mask, ...
    'bwmask_ca', mask_ca, ...
    'boundary', {{}}, ...
    'boundary_ca', {{}}, ...
    'position', {{}}, ...
    'position_ca', {{}});
traces = [];
traces_ca = [];

nrois = max(max(mask(:)), max(mask_ca(:)));
for roi_idx = 1:nrois
    mask_v = mask == roi_idx;
    mask_c = mask_ca == roi_idx;
    if ~any(mask_v(:)) || ~any(mask_c(:))
        continue;
    end

    boundary_v = extract_dual_boundary(mask_v);
    boundary_ca = extract_dual_boundary(mask_c);
    color = colors(mod(roi_idx - 1, size(colors, 1)) + 1, :);

    rois.boundary{end+1} = boundary_v; %#ok<AGROW>
    rois.boundary_ca{end+1} = boundary_ca; %#ok<AGROW>
    rois.position{end+1} = [boundary_v(:, 2), boundary_v(:, 1)]; %#ok<AGROW>
    rois.position_ca{end+1} = [boundary_ca(:, 2), boundary_ca(:, 1)]; %#ok<AGROW>

    traces(:, end+1) = mean(movie(mask_v(:), :), 1, 'double')'; %#ok<AGROW>
    traces_ca(:, end+1) = mean(movie_ca(mask_c(:), :), 1, 'double')'; %#ok<AGROW>

    plot_dual_roi_overlay(axes_info, boundary_v, boundary_ca, roi_idx, color);
end

if size(mask, 1) ~= image_height || size(mask, 2) ~= image_width
    warning('Voltage ROI mask size does not match current image geometry.');
end
if size(mask_ca, 1) ~= image_height || size(mask_ca, 2) ~= image_width
    warning('Calcium ROI mask size does not match current image geometry.');
end

if ~isempty(traces)
    cla(axes_info.trace_v);
    plot(axes_info.trace_v, traces(:, end), 'Color', colors(mod(size(traces, 2) - 1, size(colors, 1)) + 1, :), 'LineWidth', 1);
    xlabel(axes_info.trace_v, 'Frames');
    ylabel(axes_info.trace_v, 'Voltage Intensity');
    grid(axes_info.trace_v, 'on');
    title(axes_info.trace_v, sprintf('Voltage ROI %d', size(traces, 2)));

    cla(axes_info.trace_ca);
    plot(axes_info.trace_ca, traces_ca(:, end), 'Color', colors(mod(size(traces_ca, 2) - 1, size(colors, 1)) + 1, :), 'LineWidth', 1);
    xlabel(axes_info.trace_ca, 'Frames');
    ylabel(axes_info.trace_ca, 'Calcium Intensity');
    grid(axes_info.trace_ca, 'on');
    title(axes_info.trace_ca, sprintf('Calcium ROI %d', size(traces_ca, 2)));
end
end

function boundary = extract_dual_boundary(mask)
boundary_cells = bwboundaries(mask);
if isempty(boundary_cells)
    boundary = zeros(0, 2);
else
    boundary = boundary_cells{1};
end
end

function [mask_voltage, mask_calcium, note_text] = reconcile_reused_dual_masks(mask_voltage, mask_calcium, x_offset, y_offset, offset_mode)
note_text = "";
if isempty(mask_voltage) && isempty(mask_calcium)
    return;
end

rebuild_from_voltage = ~isempty(mask_voltage) && ~strcmpi(offset_mode, "none");
rebuild_from_calcium = isempty(mask_voltage) && ~isempty(mask_calcium) && ~strcmpi(offset_mode, "none");

if rebuild_from_voltage
    mask_calcium = translate_dual_label_mask(mask_voltage, -x_offset, -y_offset);
    note_text = sprintf(['Reuse ROI path: calcium mask was rebuilt from the saved voltage ROI mask ' ...
        'using offset [%.3f %.3f].'], x_offset, y_offset);
    return;
end

if rebuild_from_calcium
    mask_voltage = translate_dual_label_mask(mask_calcium, x_offset, y_offset);
    note_text = sprintf(['Reuse ROI path: voltage mask was rebuilt from the saved calcium ROI mask ' ...
        'using offset [%.3f %.3f].'], x_offset, y_offset);
    return;
end

if isempty(mask_voltage) && ~isempty(mask_calcium)
    note_text = 'Reuse ROI path: only calcium mask is available and offset mode is none.';
elseif ~isempty(mask_voltage) && isempty(mask_calcium)
    note_text = 'Reuse ROI path: only voltage mask is available and offset mode is none.';
end
end

function translated_mask = translate_dual_label_mask(mask_in, dx, dy)
mask_in = double(mask_in);
if exist('imtranslate', 'file') == 2
    translated_mask = imtranslate(mask_in, [dx, dy], 'nearest', 'FillValues', 0);
else
    translated_mask = shift_dual_label_mask_integer(mask_in, dx, dy);
end
translated_mask = round(translated_mask);
translated_mask(translated_mask < 0) = 0;
end

function shifted_mask = shift_dual_label_mask_integer(mask_in, dx, dy)
shift_x = round(dx);
shift_y = round(dy);
shifted_mask = zeros(size(mask_in));
[nrows_mask, ncols_mask] = size(mask_in);

src_rows = max(1, 1 - shift_y):min(nrows_mask, nrows_mask - shift_y);
src_cols = max(1, 1 - shift_x):min(ncols_mask, ncols_mask - shift_x);
dst_rows = src_rows + shift_y;
dst_cols = src_cols + shift_x;
if isempty(src_rows) || isempty(src_cols)
    return;
end
shifted_mask(dst_rows, dst_cols) = mask_in(src_rows, src_cols);
end

function [rois, traces, traces_ca, overlay_handles] = remove_last_dual_roi(rois, traces, traces_ca, overlay_handles, axes_info)
last_roi = max(rois.bwmask(:));
if last_roi <= 0
    return;
end

rois.bwmask(rois.bwmask == last_roi) = 0;
rois.bwmask_ca(rois.bwmask_ca == last_roi) = 0;
if ~isempty(rois.boundary)
    rois.boundary(end) = [];
end
if ~isempty(rois.boundary_ca)
    rois.boundary_ca(end) = [];
end
if ~isempty(rois.position)
    rois.position(end) = [];
end
if ~isempty(rois.position_ca)
    rois.position_ca(end) = [];
end
if ~isempty(traces)
    traces(:, end) = [];
end
if ~isempty(traces_ca)
    traces_ca(:, end) = [];
end

if ~isempty(overlay_handles)
    handles = overlay_handles{end};
    handles = handles(isgraphics(handles));
    delete(handles);
    overlay_handles(end) = [];
end

cla(axes_info.trace_v);
grid(axes_info.trace_v, 'on');
xlabel(axes_info.trace_v, 'Frames');
ylabel(axes_info.trace_v, 'Voltage Intensity');
if ~isempty(traces)
    plot(axes_info.trace_v, traces(:, end), 'k', 'LineWidth', 1);
    title(axes_info.trace_v, sprintf('Voltage ROI %d', size(traces, 2)));
else
    title(axes_info.trace_v, 'Voltage Trace');
end

cla(axes_info.trace_ca);
grid(axes_info.trace_ca, 'on');
xlabel(axes_info.trace_ca, 'Frames');
ylabel(axes_info.trace_ca, 'Calcium Intensity');
if ~isempty(traces_ca)
    plot(axes_info.trace_ca, traces_ca(:, end), 'k', 'LineWidth', 1);
    title(axes_info.trace_ca, sprintf('Calcium ROI %d', size(traces_ca, 2)));
else
    title(axes_info.trace_ca, 'Calcium Trace');
end
end

function [x_offset, y_offset] = estimate_dual_offset(axes_info, mean_voltage, mean_calcium, image_height, image_width, offset_mode)
switch lower(string(offset_mode))
    case "manual_points"
        [x_offset, y_offset] = estimate_dual_offset_from_manual_points(axes_info);
    case "matlab_register"
        [x_offset, y_offset] = estimate_dual_offset_with_matlab_register( ...
            axes_info, mean_voltage, mean_calcium, image_height, image_width);
    otherwise
        error('Unsupported offset_mode inside select_ROI_dual: %s', char(string(offset_mode)));
end
fprintf('Dual ROI offset estimated (%s): [x y] = [%.3f %.3f]\n', ...
    char(string(offset_mode)), x_offset, y_offset);
end

function [x_offset, y_offset] = estimate_dual_offset_from_manual_points(axes_info)
disp(['Manual dual offset estimation: click one point on the Voltage Image ' ...
    'and the matching point on the Calcium Image. ' ...
    'The resulting offset moves calcium into voltage coordinates.']);

voltage_point = request_dual_point(axes_info.image_v, [0.85, 0.2, 0.2], 'Voltage point');
calcium_point = request_dual_point(axes_info.image_ca, [0.2, 0.7, 0.2], 'Calcium point');

x_offset = voltage_point(1) - calcium_point(1);
y_offset = voltage_point(2) - calcium_point(2);
end

function [x_offset, y_offset] = estimate_dual_offset_with_matlab_register(axes_info, mean_voltage, mean_calcium, image_height, image_width)
if exist('imregconfig', 'file') ~= 2 || exist('imregtform', 'file') ~= 2
    error(['MATLAB built-in registration requires imregconfig and imregtform ' ...
        '(Image Processing Toolbox).']);
end

fixed_image = normalize_dual_registration_image(mean_voltage, image_height, image_width);
moving_image = normalize_dual_registration_image(mean_calcium, image_height, image_width);

% Registration direction is fixed on purpose:
%   moving = calcium
%   fixed  = voltage
% Therefore the returned translation moves calcium into voltage space,
% matching the offset convention used by ROI projection.
[optimizer, metric] = imregconfig('multimodal');
tform = imregtform(moving_image, fixed_image, 'translation', optimizer, metric);
x_offset = double(tform.T(3, 1));
y_offset = double(tform.T(3, 2));

show_dual_registration_preview(axes_info, fixed_image, moving_image, [x_offset, y_offset]);
end

function point_position = request_dual_point(ax, color, point_label)
axes(ax); %#ok<LAXES>
point_handle = drawpoint('Parent', ax, 'Color', color);
point_position = double(point_handle.Position);
if any(~isfinite(point_position)) || numel(point_position) ~= 2
    delete_if_graphics(point_handle);
    error('Failed to capture %s during manual offset estimation.', point_label);
end
label_handle = text(ax, point_position(1), point_position(2), [' ', point_label], ...
    'Color', color, 'FontWeight', 'bold', 'Interpreter', 'none');
drawnow;
pause(0.05);
delete_if_graphics([point_handle, label_handle]);
end

function image_out = normalize_dual_registration_image(image_in, image_height, image_width)
image_out = double(image_in);
if ~isequal(size(image_out), [image_height, image_width])
    error('Registration image size mismatch: expected [%d %d].', image_height, image_width);
end
image_out(~isfinite(image_out)) = 0;
if exist('mat2gray', 'file') == 2
    image_out = mat2gray(image_out);
else
    finite_pixels = image_out(isfinite(image_out));
    if isempty(finite_pixels)
        image_out = zeros(size(image_out));
    else
        pixel_min = min(finite_pixels);
        pixel_max = max(finite_pixels);
        if pixel_max > pixel_min
            image_out = (image_out - pixel_min) ./ (pixel_max - pixel_min);
        else
            image_out = zeros(size(image_out));
        end
    end
end
end

function show_dual_registration_preview(axes_info, fixed_image, moving_image, offset_xy)
registered_moving = moving_image;
if exist('imtranslate', 'file') == 2
    registered_moving = imtranslate(moving_image, offset_xy, 'FillValues', 0);
end

cla(axes_info.image_v);
imshowpair(fixed_image, registered_moving, 'Parent', axes_info.image_v);
title(axes_info.image_v, sprintf('Voltage vs Registered Calcium | [%.2f %.2f]', offset_xy(1), offset_xy(2)));

cla(axes_info.image_ca);
imshowpair(fixed_image, moving_image, 'Parent', axes_info.image_ca);
title(axes_info.image_ca, 'Voltage vs Raw Calcium');
end

function offset_mode = normalize_dual_offset_mode(offset_mode)
if islogical(offset_mode) || (isnumeric(offset_mode) && isscalar(offset_mode))
    if logical(offset_mode)
        offset_mode = "manual_points";
    else
        offset_mode = "none";
    end
    return;
end

offset_mode = lower(strtrim(string(offset_mode)));
switch offset_mode
    case {"none", "off", "false", "0", "reuse_or_zero"}
        offset_mode = "none";
    case {"manual", "manual_point", "manual_points", "points", "point"}
        offset_mode = "manual_points";
    case {"matlab_register", "register", "imreg", "imregtform"}
        offset_mode = "matlab_register";
    otherwise
        error(['Unsupported dual offset mode: %s. Use ''none'', ' ...
            '''manual_points'', or ''matlab_register''.'], char(offset_mode));
end
end

function set_dual_key_pressed(event, fig)
if any(strcmp(event.Key, {'space', 'return', 'r', 'q', 'o'}))
    fig.UserData.dual_key = event.Key;
end
end

function key = wait_for_dual_key(fig)
fig.UserData.dual_key = [];
waitfor(fig, 'UserData');
if isfield(fig.UserData, 'dual_key')
    key = fig.UserData.dual_key;
else
    key = '';
end
end

function delete_if_graphics(handle)
if ~isempty(handle) && all(isgraphics(handle))
    delete(handle);
end
end

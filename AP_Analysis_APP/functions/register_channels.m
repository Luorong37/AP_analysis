function [registration, preview] = register_channels( ...
        voltage_movie, calcium_movie, voltage_profile, calcium_profile, spec)
%REGISTER_CHANNELS Estimate the calcium-to-voltage coordinate translation.
% OFFSET_XY follows the frozen Dual convention:
%   voltage_position = calcium_position + [x y].
% Movies are never moved by this function. PREVIEW contains normalized mean
% images for UI or file rendering by the caller.
%
% SPEC fields:
%   mode          "none", "manual_points", or "matlab_register"
%   reuse_offset  optional [x y] used by mode "none"
%   manual_points optional [voltage_x voltage_y; calcium_x calcium_y]
% Input formats:
%   voltage_movie/calcium_movie Numeric/logical movies accepted by
%       normalize_movie_matrix; analysis frame sizes must match exactly.
%   voltage_profile/calcium_profile Bound profiles carrying their roles and
%       roi.frame_size; neither profile contains a request.
% Outputs:
%   registration Scalar struct with offset_xy [x y] in pixel coordinates.
%   preview Scalar struct of normalized [ncols x nrows] mean images.

if nargin < 5 || isempty(spec)
    spec = struct('mode',"none",'reuse_offset',[]);
end
validate_spec(spec);

[voltage_matrix, voltage_info] = normalize_movie_matrix( ...
    voltage_movie, voltage_profile);
[calcium_matrix, calcium_info] = normalize_movie_matrix( ...
    calcium_movie, calcium_profile);
if ~isequal(voltage_info.frame_size, calcium_info.frame_size)
    error('AAA:Functions:RegistrationGeometryMismatch', ...
        'Voltage and calcium analysis frame sizes must match exactly.');
end

frame_size = voltage_info.frame_size;
fixed = normalize_image(reshape(mean(double(voltage_matrix), 2), frame_size));
moving = normalize_image(reshape(mean(double(calcium_matrix), 2), frame_size));
preview = struct('voltage_mean',fixed,'calcium_mean',moving, ...
    'calcium_registered',moving);

mode = lower(string(spec.mode));
offset_xy = [];
applied = false;
note = "";
switch mode
    case "none"
        if isfield(spec, 'reuse_offset') && ~isempty(spec.reuse_offset)
            offset_xy = double(reshape(spec.reuse_offset, 1, 2));
            applied = any(offset_xy ~= 0);
            note = "Used the supplied saved offset; no new registration was estimated.";
        else
            offset_xy = [0 0];
            note = "No offset was estimated or supplied; [0 0] is used.";
        end
    case "manual_points"
        if isfield(spec, 'manual_points') && ~isempty(spec.manual_points)
            points = double(spec.manual_points);
        else
            points = pick_manual_points(fixed, moving);
        end
        offset_xy = points(1,:) - points(2,:);
        applied = true;
        note = "Offset was picked from matching voltage and calcium points.";
    case "matlab_register"
        if exist('imregconfig', 'file') ~= 2 || exist('imregtform', 'file') ~= 2
            error('AAA:Functions:RegistrationToolboxMissing', ...
                ['matlab_register requires imregconfig and imregtform ' ...
                 'from Image Processing Toolbox.']);
        end
        [optimizer, metric] = imregconfig('multimodal');
        tform = imregtform(moving, fixed, 'translation', optimizer, metric);
        offset_xy = [double(tform.T(3,1)), double(tform.T(3,2))];
        applied = true;
        note = "Estimated a translation mapping calcium into voltage coordinates.";
end

preview.calcium_registered = translate_image(moving, offset_xy);
registration = struct( ...
    'mode', mode, ...
    'applied', applied, ...
    'offset_xy', offset_xy, ...
    'coordinate_rule', "voltage_position = calcium_position + offset_xy", ...
    'movie_transform', "none", ...
    'note', note, ...
    'created_at', datetime("now"));
end

function validate_spec(spec)
if ~isstruct(spec) || ~isscalar(spec) || ~isfield(spec, 'mode')
    error('AAA:Functions:InvalidRegistrationSpec', ...
        'Registration spec must be one scalar struct containing mode.');
end
mode = lower(string(spec.mode));
if ~isscalar(mode) || ~ismember(mode, ...
        ["none","manual_points","matlab_register"])
    error('AAA:Functions:InvalidRegistrationMode', ...
        'Registration mode must be none, manual_points, or matlab_register.');
end
if isfield(spec, 'reuse_offset') && ~isempty(spec.reuse_offset)
    validate_offset(spec.reuse_offset, 'reuse_offset');
end
if isfield(spec, 'manual_points') && ~isempty(spec.manual_points)
    value = spec.manual_points;
    if ~isnumeric(value) || ~isequal(size(value), [2 2]) ...
            || any(~isfinite(value), 'all')
        error('AAA:Functions:InvalidManualPoints', ...
            'manual_points must be a finite 2-by-2 [x y] matrix.');
    end
end
end

function validate_offset(value, label)
if ~isnumeric(value) || numel(value) ~= 2 || any(~isfinite(value))
    error('AAA:Functions:InvalidRegistrationOffset', ...
        '%s must contain two finite numeric [x y] values.', label);
end
end

function points = pick_manual_points(fixed, moving)
if ~usejava('desktop')
    error('AAA:Functions:ManualRegistrationRequiresDesktop', ...
        'manual_points without supplied points requires MATLAB Desktop.');
end
fig = figure('Color','w','Name','AAA Manual Channel Registration', ...
    'Position',[100 100 1200 520]);
cleanup = onCleanup(@() close_if_valid(fig)); %#ok<NASGU>
ax1 = subplot(1,2,1,'Parent',fig);
imshow(fixed,'Parent',ax1);
title(ax1,'Voltage average: click reference point');
ax2 = subplot(1,2,2,'Parent',fig);
imshow(moving,'Parent',ax2);
title(ax2,'Calcium average: click matching point');
axes(ax1); %#ok<LAXES>
[x1,y1] = ginput(1);
axes(ax2); %#ok<LAXES>
[x2,y2] = ginput(1);
points = [x1 y1; x2 y2];
end

function image_out = normalize_image(image_in)
image_out = double(image_in);
image_out(~isfinite(image_out)) = 0;
minimum = min(image_out, [], 'all');
maximum = max(image_out, [], 'all');
if maximum > minimum
    image_out = (image_out - minimum) ./ (maximum - minimum);
else
    image_out = zeros(size(image_out));
end
end

function shifted = translate_image(image_in, offset_xy)
validate_offset(offset_xy, 'offset_xy');
if exist('imtranslate', 'file') == 2
    shifted = imtranslate(double(image_in), offset_xy, 'FillValues', 0);
    return;
end
shifted = zeros(size(image_in));
dx = round(offset_xy(1));
dy = round(offset_xy(2));
rows = max(1,1-dy):min(size(image_in,1),size(image_in,1)-dy);
cols = max(1,1-dx):min(size(image_in,2),size(image_in,2)-dx);
if ~isempty(rows) && ~isempty(cols)
    shifted(rows+dy,cols+dx) = image_in(rows,cols);
end
end

function close_if_valid(fig)
if isgraphics(fig)
    close(fig);
end
end

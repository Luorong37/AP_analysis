function [quick_map, info] = create_map(movie, profile)
%CREATE_MAP Create the shared Analysis3 voltage/calcium activity map.
% PROFILE is the complete bound channel profile. This function performs no
% plotting, file I/O, ROI selection, truncation, padding, or interpolation.
% Inputs:
%   movie Numeric/logical [ncols x nrows x frames], one spatial frame, or
%         [pixels x frames] using MATLAB linear pixel order.
%   profile Scalar bound profile requiring role, roi.frame_size
%           [ncols nrows], and positive integer roi.map.bin.
% Outputs:
%   quick_map Double [ncols x nrows] activity map.
%   info Scalar struct recording role, geometry, and frozen map rules.

validate_profile(profile);
[movie_matrix, movie_info] = normalize_movie_matrix(movie, profile);
frame_size = movie_info.frame_size;
map_bin = double(profile.roi.map.bin);
if any(mod(frame_size, map_bin) ~= 0)
    error('AAA:Functions:MapBinGeometryMismatch', ...
        ['profile.roi.frame_size [%d %d] must be divisible by map bin %d. ' ...
         'The map function does not crop, pad, or reshape across frames.'], ...
        frame_size(1), frame_size(2), map_bin);
end

nframes = movie_info.frame_count;
movie_3d = reshape(movie_matrix, frame_size(1), frame_size(2), nframes);
kernel = ones(3,3);
kernel(5) = 0;
for frame_idx = 1:nframes
    movie_3d(:,:,frame_idx) = imfilter( ...
        movie_3d(:,:,frame_idx) ./ 8, kernel, 'conv');
end

filtered = reshape(movie_3d, prod(frame_size), nframes);
quick_map = zeros(prod(frame_size), 1);
role = lower(string(profile.role));
for pixel_idx = 1:size(filtered, 1)
    switch role
        case "voltage"
            corrected = detrend(double(filtered(pixel_idx,:)), 1);
            maximum = max(abs(corrected));
            index = find(abs(corrected) == maximum, 1, 'first');
            quick_map(pixel_idx) = corrected(index);
        case "calcium"
            corrected = detrend(double(filtered(pixel_idx,:)), 2);
            quick_map(pixel_idx) = std(corrected);
    end
end

quick_map = reshape(quick_map, frame_size);
if any(quick_map > 0, 'all')
    quick_map(quick_map > 0) = quick_map(quick_map > 0) ...
        - mean(quick_map(quick_map > 0));
end
% Preserve the frozen Dual order exactly: the negative mask is evaluated
% after positive values have been centered, rather than cached beforehand.
if any(quick_map < 0, 'all')
    quick_map(quick_map < 0) = quick_map(quick_map < 0) ...
        - mean(quick_map(quick_map < 0));
end

info = struct( ...
    'role', role, ...
    'frame_size', frame_size, ...
    'frame_count', nframes, ...
    'map_bin', map_bin, ...
    'spatial_filter', "3x3 neighbors excluding center, divided by 8", ...
    'voltage_rule', "signed largest absolute first-order detrended residual", ...
    'calcium_rule', "standard deviation of second-order detrended trace", ...
    'centering_rule', "positive and negative pixels centered separately", ...
    'geometry_rule', "exact divisibility; no crop, pad, or interpolation");
end

function validate_profile(profile)
if ~isstruct(profile) || ~isscalar(profile) ...
        || ~all(isfield(profile, {'role','roi'})) ...
        || ~isstruct(profile.roi) || ~isscalar(profile.roi) ...
        || ~all(isfield(profile.roi, {'frame_size','map'})) ...
        || ~isstruct(profile.roi.map) || ~isscalar(profile.roi.map) ...
        || ~isfield(profile.roi.map, 'bin')
    error('AAA:Functions:InvalidMapProfile', ...
        'Map requires profile.role, profile.roi.frame_size, and profile.roi.map.bin.');
end
role = lower(string(profile.role));
if ~isscalar(role) || ~ismember(role, ["voltage","calcium"])
    error('AAA:Functions:InvalidMapRole', ...
        'profile.role must be voltage or calcium.');
end
bin = profile.roi.map.bin;
if ~isnumeric(bin) || ~isscalar(bin) || ~isfinite(bin) ...
        || bin < 1 || bin ~= round(bin)
    error('AAA:Functions:InvalidMapBin', ...
        'profile.roi.map.bin must be a positive integer.');
end
end

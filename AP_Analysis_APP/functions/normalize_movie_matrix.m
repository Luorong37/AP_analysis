function [movie_matrix, info] = normalize_movie_matrix(movie, profile)
%NORMALIZE_MOVIE_MATRIX Validate one channel movie without changing values.
% PROFILE is the complete bound channel profile. The role has one authority
% at profile.role; runtime geometry is stored at profile.roi.frame_size.
% No spatial or temporal truncation, padding, interpolation, or casting is
% performed. The returned matrix is [pixels x frames] in MATLAB linear order.
% Inputs:
%   movie Numeric/logical [ncols x nrows x frames], one
%         [ncols x nrows] spatial frame, or [pixels x frames].
%   profile Scalar bound profile requiring role and roi.frame_size
%           [ncols nrows]. It contains no request.
% Outputs:
%   movie_matrix Same-class [pixels x frames] matrix.
%   info Scalar struct describing the recognized input layout and geometry.

require_profile(profile);
role = lower(string(profile.role));
frame_size = double(reshape(profile.roi.frame_size, 1, 2));
expected_pixels = prod(frame_size);

if (~isnumeric(movie) && ~islogical(movie)) || isempty(movie) ...
        || ~isreal(movie) || ndims(movie) > 3
    error('AAA:Algorithms:InvalidRoiMovie', ...
        'ROI algorithms require one nonempty real 2-D or 3-D numeric movie.');
end

if ndims(movie) == 3
    if ~isequal([size(movie, 1), size(movie, 2)], frame_size)
        error('AAA:Algorithms:RoiMovieGeometryMismatch', ...
            'Movie geometry [%d %d] does not match profile.roi.frame_size [%d %d].', ...
            size(movie, 1), size(movie, 2), frame_size(1), frame_size(2));
    end
    nframes = size(movie, 3);
    input_layout = "spatial_3d";
elseif isequal(size(movie), frame_size)
    nframes = 1;
    input_layout = "single_spatial_frame";
elseif size(movie, 1) == expected_pixels
    nframes = size(movie, 2);
    input_layout = "pixels_by_frames";
else
    error('AAA:Algorithms:RoiMovieGeometryMismatch', ...
        ['A 2-D ROI movie must be either one [%d %d] spatial frame or ' ...
         'a [%d x frames] pixel matrix.'], ...
        frame_size(1), frame_size(2), expected_pixels);
end

movie_matrix = reshape(movie, expected_pixels, nframes);
info = struct( ...
    'role', role, ...
    'frame_size', frame_size, ...
    'frame_count', nframes, ...
    'input_layout', input_layout, ...
    'input_class', string(class(movie)), ...
    'output_class', string(class(movie_matrix)), ...
    'alignment_rule', "exact; no truncation, padding, or interpolation");
end

function require_profile(profile)
if ~isstruct(profile) || ~isscalar(profile)
    error('AAA:Algorithms:InvalidRoiProfile', ...
        'ROI algorithm profile must be one scalar struct.');
end
required = {'role', 'roi'};
missing = required(~isfield(profile, required));
if ~isempty(missing)
    error('AAA:Algorithms:InvalidRoiProfile', ...
        'ROI algorithm profile requires: %s.', strjoin(missing, ', '));
end
role = lower(string(profile.role));
if ~isscalar(role) || ~ismember(role, ["voltage", "calcium"])
    error('AAA:Algorithms:InvalidRoiRole', ...
        'profile.role must be voltage or calcium.');
end
if ~isstruct(profile.roi) || ~isscalar(profile.roi) ...
        || ~isfield(profile.roi, 'frame_size')
    error('AAA:Algorithms:InvalidRoiProfile', ...
        'profile.roi.frame_size is required.');
end
frame_size = profile.roi.frame_size;
if ~isnumeric(frame_size) || numel(frame_size) ~= 2 ...
        || any(~isfinite(frame_size)) || any(frame_size < 1) ...
        || any(frame_size ~= round(frame_size))
    error('AAA:Algorithms:InvalidRoiFrameSize', ...
        'profile.roi.frame_size must be [ncols nrows] positive integers.');
end
end

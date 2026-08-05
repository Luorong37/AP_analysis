function model = estimate_motion(movie_3d, profile)
%ESTIMATE_MOTION Estimate one rigid NoRMCorre shift field.
%   The function performs no saving and knows nothing about AP or Dual
%   workflows. Cross-channel sharing is decided by the calling section.
% Inputs:
%   movie_3d Nonempty numeric/logical [ncols x nrows x frames] movie.
%   profile Scalar bound profile requiring role and profile.motion fields
%           bin_width, max_shift, upsample_factor, iterations,
%           correct_bidir, and highpass. It contains no request.
% Output:
%   model Scalar NoRMCorre model. model.shifts has exactly one entry per
%         input frame; no temporal padding, truncation, or interpolation.

validate_movie(movie_3d);
[role, params] = resolve_profile(profile);
[ncols, nrows, nframes] = size(movie_3d);

options = NoRMCorreSetParms( ...
    'd1', ncols, ...
    'd2', nrows, ...
    'bin_width', params.bin_width, ...
    'max_shift', params.max_shift, ...
    'us_fac', params.upsample_factor, ...
    'iter', params.iterations, ...
    'correct_bidir', params.correct_bidir);

movie_single = single(movie_3d);
if params.highpass
    movie_for_estimation = create_temp_highpass(movie_single);
else
    movie_for_estimation = movie_single;
end
[~, shifts, ~] = normcorre_batch(movie_for_estimation, options);

if numel(shifts) ~= nframes
    error('AAA:Algorithms:EstimatedShiftFrameMismatch', ...
        'NoRMCorre returned %d shifts for a %d-frame movie.', ...
        numel(shifts), nframes);
end
model = struct( ...
    'method', "NoRMCorre_rigid", ...
    'role', role, ...
    'shifts', shifts, ...
    'options', options, ...
    'frame_count', nframes, ...
    'frame_size', [ncols, nrows], ...
    'parameters', params, ...
    'source_shift_file', "", ...
    'reused', false, ...
    'created_at', datetime('now'));
end

function validate_movie(movie_3d)
if ~isnumeric(movie_3d) && ~islogical(movie_3d)
    error('AAA:Algorithms:InvalidMotionMovie', ...
        'Motion input must be a numeric or logical movie.');
end
if isempty(movie_3d) || ndims(movie_3d) > 3
    error('AAA:Algorithms:InvalidMotionMovie', ...
        'Motion input must be a nonempty 2-D or 3-D movie.');
end
end

function [role, params] = resolve_profile(profile)
if ~isstruct(profile) || ~isscalar(profile) ...
        || ~isfield(profile, 'role') || strlength(string(profile.role)) == 0
    error('AAA:Algorithms:InvalidMotionProfile', ...
        'Motion profile must be scalar and contain one nonempty role.');
end
if ~isfield(profile, 'motion') || ~isstruct(profile.motion) ...
        || ~isscalar(profile.motion)
    error('AAA:Algorithms:IncompleteMotionProfile', ...
        'Bound channel profile must contain one scalar motion struct.');
end
role = string(profile.role);
params = profile.motion;
required = {'bin_width','max_shift','upsample_factor', ...
    'iterations','correct_bidir','highpass'};
for idx = 1:numel(required)
    name = required{idx};
    if ~isfield(params, name) || isempty(params.(name))
        error('AAA:Algorithms:IncompleteMotionProfile', ...
            'profile.motion is missing %s.', name);
    end
end

validate_positive_integer(params.bin_width, 'bin_width');
validate_nonnegative_integer(params.max_shift, 'max_shift');
validate_positive_integer(params.upsample_factor, 'upsample_factor');
validate_positive_integer(params.iterations, 'iterations');
if ~(islogical(params.correct_bidir) && isscalar(params.correct_bidir))
    error('AAA:Algorithms:InvalidMotionParameter', ...
        'profile.motion.correct_bidir must be one logical value.');
end
if ~(islogical(params.highpass) && isscalar(params.highpass))
    error('AAA:Algorithms:InvalidMotionParameter', ...
        'profile.motion.highpass must be one logical value.');
end
end

function validate_positive_integer(value, name)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) ...
        || value < 1 || value ~= round(value)
    error('AAA:Algorithms:InvalidMotionParameter', ...
        'profile.motion.%s must be a positive integer.', name);
end
end

function validate_nonnegative_integer(value, name)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) ...
        || value < 0 || value ~= round(value)
    error('AAA:Algorithms:InvalidMotionParameter', ...
        'profile.motion.%s must be a nonnegative integer.', name);
end
end

function movie_highpass = create_temp_highpass(movie_3d)
% Preserve the established Dual high-pass kernel used only for estimation.
gSig = 7;
gSiz = 17;
psf = fspecial('gaussian', round(2 * gSiz), gSig);
inside = psf >= max(psf(:, 1));
psf = psf - mean(psf(inside));
psf(~inside) = 0;
movie_highpass = imfilter(movie_3d, psf, 'symmetric');
end

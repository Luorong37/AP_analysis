function tests = test_aaa_channel_input_motion
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
unit_dir = fileparts(mfilename('fullpath'));
app_dir = fileparts(fileparts(unit_dir));
old_path = path;
addpath(app_dir);
AAA_startup();
testCase.TestData.oldPath = old_path;
testCase.TestData.tempRoot = tempname;
mkdir(testCase.TestData.tempRoot);
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
if isfolder(testCase.TestData.tempRoot)
    rmdir(testCase.TestData.tempRoot, 's');
end
end

function testLoadAlgorithmAppliesActualPageTranspose(testCase)
movie = reshape(uint16(1:24), [2, 3, 4]);
source_file = fullfile(testCase.TestData.tempRoot, 'transpose_source.mat');
save(source_file, 'movie');
profile = make_profile("voltage", 20, true, "");

[loaded, info] = load_movie(source_file, profile);

verifyEqual(testCase, loaded, pagetranspose(movie));
verifyEqual(testCase, info.original_frame_size, [2, 3]);
verifyEqual(testCase, info.analysis_frame_size, [3, 2]);
verifyEqual(testCase, info.frame_count, 4);
verifyTrue(testCase, info.transpose_before_analysis);
end

function testInputPreservesRolesCameraIndicesAndIndependentTimeAxes(testCase)
voltage_movie = reshape(uint16(1:18), [2, 3, 3]);
calcium_movie = reshape(uint16(1:120), [4, 6, 5]);
voltage_file = fullfile(testCase.TestData.tempRoot, 'voltage_source.mat');
calcium_file = fullfile(testCase.TestData.tempRoot, 'calcium_source.mat');
movie = voltage_movie;
save(voltage_file, 'movie');
movie = calcium_movie;
save(calcium_file, 'movie');

channels(1) = make_channel(2, voltage_file, ...
    make_profile("voltage", 30, false, ""));
channels(2) = make_channel(1, calcium_file, ...
    make_profile("calcium", 10, false, ""));
ctx = make_test_context(channels);

[ctx, outcome] = aaa.sections.common.input(ctx, "run");

verifyEqual(testCase, outcome.geometry.action, ...
    "downsample_calcium_to_voltage");
verifyEqual(testCase, size(ctx.channels(1).data.movie_3d), [2, 3, 3]);
verifyEqual(testCase, size(ctx.channels(2).data.movie_3d), [2, 3, 5]);
verifyEqual(testCase, ctx.channels(1).data.time, (1:3)' / 30);
verifyEqual(testCase, ctx.channels(2).data.time, (1:5)' / 10);
verifyEqual(testCase, ctx.channels(1).data.movie_info.camera_index, 2);
verifyEqual(testCase, ctx.channels(2).data.movie_info.camera_index, 1);
verifyEqual(testCase, ctx.channels(1).profile.role, "voltage");
verifyEqual(testCase, ctx.channels(2).profile.role, "calcium");
verifyFalse(testCase, isfield(ctx.channels, 'role'));
verifyFalse(testCase, isfield(ctx.channels, 'frame_rate'));
end

function testCreateContextResolvesCycleCameraSourcesByCameraIndex(testCase)
cycle_path = fullfile(testCase.TestData.tempRoot, 'Cycle1');
camera1_path = fullfile(cycle_path, 'Cam1_calcium');
camera2_path = fullfile(cycle_path, 'Cam2_voltage');
mkdir(camera1_path);
mkdir(camera2_path);
% Default Dual profile transposes calcium, so its raw frame is 3-by-2.
movie = reshape(uint16(1:24), [3, 2, 4]);
save(fullfile(camera1_path, 'movie.mat'), 'movie');
movie = reshape(uint16(101:124), [2, 3, 4]);
save(fullfile(camera2_path, 'movie.mat'), 'movie');
request = aaa.schema.default_request("dual", "single");
request.input_path = string(cycle_path);
ctx = aaa.sections.create_context(request, struct());

[ctx, outcome] = aaa.sections.common.input(ctx, "run");

verifyEqual(testCase, outcome.input_layout.mode, ...
    "resolved_cycle_sources");
verifyEqual(testCase, ctx.channels(1).profile.role, "voltage");
verifyEqual(testCase, ctx.channels(1).camera_index, 2);
verifyTrue(testCase, contains( ...
    ctx.channels(1).data.movie_info.source_path, "Cam2_voltage"));
verifyEqual(testCase, ctx.channels(2).profile.role, "calcium");
verifyEqual(testCase, ctx.channels(2).camera_index, 1);
verifyTrue(testCase, contains( ...
    ctx.channels(2).data.movie_info.source_path, "Cam1_calcium"));
end

function testNeutralTiffWriterRoundTrip(testCase)
movie = reshape(uint16(1:24), [3, 4, 2]);
requested_path = fullfile(testCase.TestData.tempRoot, 'roundtrip.tif');
written_files = aaa.io.write_tiff_stack(movie, requested_path);
profile = make_profile("voltage", 25, false, "");

[loaded, info] = load_movie(written_files(1), profile);

verifyEqual(testCase, numel(written_files), 1);
verifyTrue(testCase, endsWith(written_files(1), ...
    "roundtrip_stack01.tif"));
verifyEqual(testCase, loaded, movie);
verifyEqual(testCase, info.frame_count, 2);
end

function testApplyMotionRejectsShiftFrameMismatch(testCase)
movie = zeros(4, 5, 3, 'single');
options = make_options([4, 5]);
model = struct('shifts', make_zero_shifts(2), 'options', options);

verifyError(testCase, ...
    @() apply_motion(movie, model), ...
    'AAA:Algorithms:MotionShiftFrameMismatch');
end

function testEqualFrameDualReuseSharesVoltageModel(testCase)
shift_file = fullfile(testCase.TestData.tempRoot, 'shared_shift.mat');
save_shift_file(shift_file, 3, [4, 5], false);
voltage_profile = make_profile("voltage", 40, false, shift_file);
calcium_profile = make_profile("calcium", 40, false, "");
channels(1) = make_loaded_channel(2, voltage_profile, ...
    reshape(single(1:60), [4, 5, 3]));
channels(2) = make_loaded_channel(1, calcium_profile, ...
    reshape(single(61:120), [4, 5, 3]));
ctx = make_test_context(channels);

[ctx, outcome] = aaa.sections.common.motion(ctx, "reuse");

verifyTrue(testCase, outcome.shared_shift);
verifyEqual(testCase, size(ctx.channels(1).data.movie_3d, 3), 3);
verifyEqual(testCase, size(ctx.channels(2).data.movie_3d, 3), 3);
verifyEqual(testCase, ...
    ctx.channels(1).data.movie_info.motion.method, ...
    'NoRMCorre_rigid_shared_voltage_reference');
verifyEqual(testCase, ...
    ctx.channels(2).data.movie_info.motion.method, ...
    'reuse_voltage_motion_shifts');
verifyEqual(testCase, ...
    ctx.channels(1).results.motion.model.source_shift_file, ...
    string(shift_file));
verifyEqual(testCase, ...
    ctx.channels(2).results.motion.model.source_shift_file, ...
    string(shift_file));
end

function testUnequalFrameReuseRequiresIndependentCalciumModel(testCase)
shift_file = fullfile(testCase.TestData.tempRoot, 'only_reference_shift.mat');
save_shift_file(shift_file, 3, [4, 5], false);
channels(1) = make_loaded_channel(2, ...
    make_profile("voltage", 40, false, shift_file), ...
    zeros(4, 5, 3, 'single'));
channels(2) = make_loaded_channel(1, ...
    make_profile("calcium", 20, false, ""), ...
    zeros(4, 5, 5, 'single'));
ctx = make_test_context(channels);

verifyError(testCase, ...
    @() aaa.sections.common.motion(ctx, "reuse"), ...
    'AAA:Sections:MissingMotionShiftFile');
end

function testUnequalFrameReusePreservesBothFrameCounts(testCase)
voltage_shift_file = fullfile(testCase.TestData.tempRoot, 'voltage_shift.mat');
calcium_shift_file = fullfile(testCase.TestData.tempRoot, 'calcium_shift.mat');
save_shift_file(voltage_shift_file, 3, [4, 5], false);
save_shift_file(calcium_shift_file, 5, [4, 5], true);
channels(1) = make_loaded_channel(2, ...
    make_profile("voltage", 40, false, voltage_shift_file), ...
    reshape(single(1:60), [4, 5, 3]));
calcium_profile = make_profile("calcium", 20, false, calcium_shift_file);
calcium_profile.motion.saved_calcium_shift_file = calcium_shift_file;
channels(2) = make_loaded_channel(1, calcium_profile, ...
    reshape(single(1:100), [4, 5, 5]));
ctx = make_test_context(channels);

[ctx, outcome] = aaa.sections.common.motion(ctx, "reuse");

verifyFalse(testCase, outcome.shared_shift);
verifyEqual(testCase, size(ctx.channels(1).data.movie_3d, 3), 3);
verifyEqual(testCase, size(ctx.channels(2).data.movie_3d, 3), 5);
verifyEqual(testCase, ctx.channels(1).data.time, (1:3)' / 40);
verifyEqual(testCase, ctx.channels(2).data.time, (1:5)' / 20);
verifyTrue(testCase, ...
    ctx.channels(2).data.movie_info.motion.frame_count_fallback);
verifyEqual(testCase, ...
    ctx.channels(2).data.movie_info.motion.method, ...
    'NoRMCorre_rigid_independent_calcium');
end

function channel = make_channel(camera_index, source_path, profile)
channel = struct( ...
    'camera_index', camera_index, ...
    'source', struct('path',string(source_path),'label',profile.role), ...
    'profile', profile, ...
    'data', struct(), ...
    'results', struct(), ...
    'output_files', struct(), ...
    'status', "pending");
end

function ctx = make_test_context(channels)
request = aaa.schema.default_request("dual", "single");
request.input_path = string(tempdir);
ctx = aaa.sections.create_context(request, struct());
ctx.channels = channels;
ctx.input_path = "";
ctx.output_path = "";
end

function channel = make_loaded_channel(camera_index, profile, movie_3d)
channel = make_channel(camera_index, "unused", profile);
channel.data.movie_3d = movie_3d;
channel.data.movie_info = struct( ...
    'role', profile.role, ...
    'camera_index', camera_index, ...
    'frame_rate', profile.frame_rate, ...
    'frame_count', size(movie_3d, 3), ...
    'analysis_frame_size', [size(movie_3d, 1), size(movie_3d, 2)], ...
    'motion', struct('applied',false), ...
    'updated_at', datetime('now'));
channel.data.time = (1:size(movie_3d, 3))' / profile.frame_rate;
end

function profile = make_profile(role, frame_rate, transpose, saved_shift_file)
profile = struct( ...
    'role', string(role), ...
    'frame_rate', frame_rate, ...
    'transpose_before_analysis', logical(transpose), ...
    'motion', struct( ...
        'bin_width', 200, ...
        'max_shift', 30, ...
        'upsample_factor', 30, ...
        'iterations', 1, ...
        'correct_bidir', false, ...
        'highpass', true, ...
        'saved_shift_file', string(saved_shift_file), ...
        'saved_calcium_shift_file', "", ...
        'save_downsampled_tif', true, ...
        'plot_metrics', true, ...
        'downsample_factor', 40, ...
        'display_quantile_low', 0.0005, ...
        'display_quantile_high', 0.99995));
end

function options = make_options(frame_size)
options = NoRMCorreSetParms( ...
    'd1', frame_size(1), ...
    'd2', frame_size(2), ...
    'bin_width', 200, ...
    'max_shift', 30, ...
    'us_fac', 30, ...
    'iter', 1, ...
    'correct_bidir', false);
end

function shifts = make_zero_shifts(nframes)
template = struct('shifts',[0 0],'shifts_up',[0 0],'diff',0);
shifts = repmat(template, nframes, 1);
end

function save_shift_file(file_path, nframes, frame_size, is_calcium)
options = make_options(frame_size);
shifts = make_zero_shifts(nframes);
if is_calcium
    shifts_c = shifts;
    options_c = options;
    save(file_path, 'shifts_c', 'options_c');
else
    shifts_r = shifts;
    options_r = options;
    save(file_path, 'shifts_r', 'options_r');
end
end

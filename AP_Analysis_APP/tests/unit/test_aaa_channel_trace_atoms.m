function tests = test_aaa_channel_trace_atoms
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
unit_dir = fileparts(mfilename('fullpath'));
app_dir = fileparts(fileparts(unit_dir));
old_path = path;
restoredefaultpath;
addpath(app_dir);
AAA_startup();
testCase.TestData.oldPath = old_path;
testCase.TestData.appDir = app_dir;
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testProductionAtomsResolveWithoutLegacyPath(testCase)
resolved = [
    string(which('remove_bleach'))
    string(which('compute_trace_metrics'))
    string(which('smooth_trace'))
    string(which('store_trace_stage'))
    string(which('fetch_trace_stage'))
    string(which('resolve_trace_stage'))
    string(which('build_section_record'))
    ];
verifyTrue(testCase, all(startsWith(resolved, ...
    string(fullfile(testCase.TestData.appDir, 'functions')))));
verifyFalse(testCase, any(contains(resolved, 'legacy')));
end

function testStoreTraceStagePreservesParentOrderAndSchema(testCase)
movie_info = struct( ...
    'source_path', "movie.tif", ...
    'transpose_before_analysis', true, ...
    'motion', struct('applied', false));
parents = {'bleach_removed', 'noise'};
data = reshape(1:6, 3, 2);
results = store_trace_stage( ...
    struct(), 'snr', data, parents, "roi.mat", movie_info, ...
    'signal_divided_by_noise_std', struct('signal_stage', 'bleach_removed'));

stage = results.trace_results.snr;
verifyEqual(testCase, stage.data, data);
verifyEqual(testCase, stage.info.parent_results, parents);
verifyEqual(testCase, fieldnames(stage.info), { ...
    'result_name'; 'parent_results'; 'roi_file'; 'movie_source'; ...
    'transpose_before_analysis'; 'motion_applied'; 'method'; ...
    'parameters'; 'created_at'});
verifyEqual(testCase, fetch_trace_stage(results, 'snr'), data);
end

function testResolveTraceStageUsesCallerPriorityWithoutFallback(testCase)
results.trace_results.raw.data = 1;
results.trace_results.raw.info = struct();
results.trace_results.bg_removed.data = 2;
results.trace_results.bg_removed.info = struct();
[data, stage] = resolve_trace_stage( ...
    results, {'bg_removed', 'raw'});
verifyEqual(testCase, data, 2);
verifyEqual(testCase, stage, 'bg_removed');

verifyError(testCase, ...
    @() resolve_trace_stage(results, {'snr', 'sensitivity'}), ...
    '');
end

function testLinearBleachMatchesFrozenFormula(testCase)
trace = [1 4; 2 3; 5 2; 9 1];
[actual, baseline, parameters] = remove_bleach( ...
    trace, struct('mode', "linear"));
expected = detrend(trace, 1);
verifyEqual(testCase, actual, expected);
verifyEqual(testCase, baseline, trace - expected);
verifyEqual(testCase, parameters, struct());
end

function testHighpassBleachMatchesFrozenFormula(testCase)
assumeEqual(testCase, exist('highpass', 'file'), 2);
frame_rate = 400;
nframes = 1000;
t_axis = (1:nframes)' / frame_rate;
trace = sin((1:nframes)' / 11) + (1:nframes)' / nframes;
spec = struct('mode', "highpass", ...
    'frame_rate', frame_rate, 'time_axis', t_axis);
[actual, baseline, parameters] = remove_bleach(trace, spec);

fc = 0.5 / t_axis(end);
padlength = round(0.05 * size(trace, 1));
padded = [repmat(trace(1, :), padlength, 1); trace; ...
    repmat(trace(end, :), padlength, 1)];
expected_padded = highpass(padded, fc, frame_rate);
expected = expected_padded(padlength+1:end-padlength, :);
verifyEqual(testCase, actual, expected);
verifyEqual(testCase, baseline, smoothdata(trace - expected));
verifyEqual(testCase, parameters, struct('fc', fc));
end

function testVoltageMetricsMatchFrozenFormulaAndEpsClamp(testCase)
assumeEqual(testCase, exist('wdenoise', 'file'), 2);
nframes = 8192;
sample = (1:nframes)';
trace = [sin(sample / 5) + 0.2 * sin(sample / 31), ...
    cos(sample / 7) - 0.1 * sin(sample / 43)];
baseline = [zeros(size(trace, 1), 1), 2 * ones(size(trace, 1), 1)];
profile = struct( ...
    'role', "voltage", ...
    'trace', struct('noise_reference', struct( ...
        'method', "wdenoise", ...
        'level', 8, ...
        'denoising_method', "FDR", ...
        'wavelet_name', "bior6.8")));
[actual, info] = compute_trace_metrics(trace, baseline, profile);

expected_ref = wdenoise(double(trace), 8, ...
    DenoisingMethod='FDR', Wavelet='bior6.8');
expected_noise = double(trace) - expected_ref;
baseline_safe = double(baseline);
baseline_safe(abs(baseline_safe) < eps) = eps;
expected_sensitivity = double(trace) ./ baseline_safe;
noise_std = std(expected_noise, 0, 1);
noise_std(noise_std < eps) = eps;
expected_snr = double(trace) ./ noise_std;

verifyEqual(testCase, actual.noise_reference, expected_ref);
verifyEqual(testCase, actual.noise, expected_noise);
verifyEqual(testCase, actual.sensitivity, expected_sensitivity);
verifyEqual(testCase, actual.snr, expected_snr);
verifyEqual(testCase, info.noise_reference_parameters, ...
    struct('level', 8, 'denoising_method', 'FDR', ...
    'wavelet_name', 'bior6.8'));
end

function testCalciumMetricsMatchFrozenFormula(testCase)
assumeEqual(testCase, exist('butter', 'file'), 2);
assumeEqual(testCase, exist('filtfilt', 'file'), 2);
nframes = 120;
trace = [(1:nframes)' + sin((1:nframes)' / 3), ...
    2 * (1:nframes)' + cos((1:nframes)' / 7)];
baseline = 100 + zeros(size(trace));
normalized_cutoff = min(0.49, 20 / max(400, 1));
profile = struct( ...
    'role', "calcium", ...
    'trace', struct('noise_reference', struct( ...
        'method', "butterworth_filtfilt", ...
        'order', 2, ...
        'normalized_cutoff', normalized_cutoff)));
[actual, info] = compute_trace_metrics(trace, baseline, profile);

[b, a] = butter(2, normalized_cutoff);
expected_ref = zeros(size(double(trace)));
for idx = 1:size(trace, 2)
    expected_ref(:, idx) = filtfilt(b, a, double(trace(:, idx)));
end
expected_noise = double(trace) - expected_ref;
expected_sensitivity = double(trace) ./ double(baseline);
noise_std = std(expected_noise, 0, 1);
noise_std(noise_std < eps) = eps;
expected_snr = double(trace) ./ noise_std;

verifyEqual(testCase, actual.noise_reference, expected_ref);
verifyEqual(testCase, actual.noise, expected_noise);
verifyEqual(testCase, actual.sensitivity, expected_sensitivity);
verifyEqual(testCase, actual.snr, expected_snr);
verifyEqual(testCase, info.noise_reference_parameters, ...
    struct('order', 2, 'normalized_cutoff', normalized_cutoff));
end

function testSmoothTraceMatchesThreeFrozenCalciumCalls(testCase)
input = reshape(1:30, 10, 3);
spec = struct('method', "movmean", 'window', 4, 'dimension', 1);
[actual, info] = smooth_trace(input, spec);
verifyEqual(testCase, actual, movmean(double(input), 4, 1));
verifyEqual(testCase, info.parameters, struct('window', 4, 'dimension', 1));
end

function testRoleAndMethodCannotDriftApart(testCase)
profile = struct( ...
    'role', "voltage", ...
    'trace', struct('noise_reference', struct( ...
        'method', "butterworth_filtfilt", ...
        'order', 2, ...
        'normalized_cutoff', 0.05)));
verifyError(testCase, ...
    @() compute_trace_metrics(ones(20, 1), ones(20, 1), profile), ...
    'AAA:Algorithms:RoleMethodMismatch');
end

function testSectionRecordMatchesFrozenSchema(testCase)
record = build_section_record( ...
    'Trace', 'description', struct('x', 1), struct('p', 2), ...
    struct('y', 3), struct('movie', 'not saved'), 'rerun note');
verifyEqual(testCase, fieldnames(record), { ...
    'section_name'; 'description'; 'input'; 'parameters'; 'output'; ...
    'excluded_inputs'; 'rerun'; 'created_at'});
verifyEqual(testCase, record.rerun.movie_saved, false);
verifyEqual(testCase, record.rerun.note, "rerun note");
end

function tests = test_aaa_peak_atoms
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

function testDualFormulaAndTableLayout(testCase)
t = (1:120)';
trace = [sin(2*pi*t/20), -sin(2*pi*t/30)];
profile = peak_profile("auto", 0.2, "absolute", 8, 0);
detected = detect_peaks(trace, profile);

verifyEqual(testCase, detected.peak_table.Properties.VariableNames, ...
    {'peak_id','roi','index','time_s','polarity', ...
    'amplitude_sensitivity','status','source','parent_peak_id', ...
    'created_stage','created_at'});
for roi_idx = 1:2
    polarity = detected.polarity{roi_idx};
    [~, expected_index] = findpeaks(trace(:,roi_idx) * polarity, ...
        'MinPeakProminence', 0.2, 'MinPeakDistance', 8, ...
        'MinPeakHeight', 0);
    verifyEqual(testCase, detected.index{roi_idx}, expected_index);
    verifyEqual(testCase, detected.amplitude{roi_idx}, ...
        trace(expected_index, roi_idx));
end
end

function testRelativeProminenceMatchesFrozenRule(testCase)
trace = [0; 1; 0; 3; 0; 2; 0];
profile = peak_profile("positive", 0.5, "relative_factor", 1, 0);
actual = detect_peaks(trace, profile);
scale = max(trace, [], 'omitnan') - mean(trace, 'omitnan');
[~, expected] = findpeaks(trace, ...
    'MinPeakProminence', 0.5 * scale, ...
    'MinPeakDistance', 1, 'MinPeakHeight', 0);
verifyEqual(testCase, actual.index{1}, expected);
end

function testRoleAndPeakParametersStayBound(testCase)
profile = peak_profile("positive", 0.1, "absolute", 1, 0);
profile.role = "calcium";
verifyError(testCase, ...
    @() detect_peaks([0;1;0], profile), ...
    'AAA:Algorithms:UnsupportedPeakRole');
end

function testTableConversionSortsAcceptedAndPreservesFallback(testCase)
profile = peak_profile("positive", 0, "absolute", 1, 0);
detected = detect_peaks([0;2;0;1;0], profile);
table_value = detected.peak_table([2 1], :);
table_value.status(1) = "deleted";
[index, amplitude, polarity] = peak_table_to_cells( ...
    table_value, 2, {1,-1});
verifyEqual(testCase, index{1}, table_value.index(2));
verifyEqual(testCase, amplitude{1}, table_value.amplitude_sensitivity(2));
verifyEqual(testCase, polarity{2}, -1);
end

function testPeakReuseValidation(testCase)
profile = peak_profile("positive", 0.1, "absolute", 1, 0);
detected = detect_peaks([0;1;0;2;0], profile);
peak_results.accepted_for_events.data = detected;
[valid, reason] = validate_peak_result(peak_results, 1, 5);
verifyTrue(testCase, valid);
verifyTrue(testCase, contains(reason, "match"));

peak_results.accepted_for_events.data.index{1}(1) = 9;
[valid, reason] = validate_peak_result(peak_results, 1, 5);
verifyFalse(testCase, valid);
verifyTrue(testCase, contains(reason, "ROI 1"));
end

function testContextBindsModeSpecificPeakDefaults(testCase)
ap_request = aaa.schema.default_request("ap", "single");
ap_request.input_path = testCase.TestData.appDir;
ap_context = aaa.sections.create_context(ap_request, struct());
verifyEqual(testCase, ap_context.channels.profile.role, "voltage");
verifyEqual(testCase, ...
    ap_context.channels.profile.peak.min_peak_distance_frames, ...
    max(2, round(0.003 * ap_request.params.ap.freq)));
verifyEqual(testCase, ...
    ap_context.channels.profile.peak.min_peak_prominence, ...
    ap_request.params.ap.peak_min_prominence);

dual_request = aaa.schema.default_request("dual", "single");
dual_request.input_path = testCase.TestData.appDir;
dual_context = aaa.sections.create_context(dual_request, struct());
profiles = [dual_context.channels.profile];
voltage_idx = find(string({profiles.role}) == "voltage", 1);
calcium_idx = find(string({profiles.role}) == "calcium", 1);
verifyEqual(testCase, ...
    profiles(voltage_idx).peak.min_peak_distance_frames, ...
    max(2, round(0.01 * profiles(voltage_idx).frame_rate)));
verifyEqual(testCase, profiles(voltage_idx).peak.global_polarity, ...
    dual_request.params.dual.voltage_polarity);
verifyTrue(testCase, isempty(fieldnames(profiles(calcium_idx).peak)));
end

function testEditorDisabledPreservesDetectedTable(testCase)
trace = [0;1;0;2;0];
profile = peak_profile("positive", 0.1, "absolute", 1, 0);
profile.peak.run_manual_edit = false;
detected = detect_peaks(trace, profile);
edited = edit_peaks(trace, detected, profile);
verifyEqual(testCase, edited.peak_table, detected.peak_table);
verifyEmpty(testCase, edited.edit_history);
verifyEqual(testCase, edited.index, detected.index);
verifyEqual(testCase, edited.amplitude, detected.amplitude);
verifyEqual(testCase, edited.polarity, detected.polarity);
end

function profile = peak_profile(mode, prominence, prominence_mode, distance, height)
profile = struct( ...
    'role', "voltage", ...
    'frame_rate', 100, ...
    'peak', struct( ...
        'polarity_mode', mode, ...
        'global_polarity', -1, ...
        'roi_polarity', [], ...
        'min_peak_prominence', prominence, ...
        'min_peak_prominence_mode', prominence_mode, ...
        'min_peak_distance_frames', distance, ...
        'min_peak_height', height, ...
        'run_manual_edit', false));
end

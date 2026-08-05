function tests = test_aaa_roi_atoms
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

function testAtomsResolveOutsideLegacy(testCase)
resolved = [
    string(which('normalize_movie_matrix'))
    string(which('extract_roi_traces'))
    string(which('remove_background'))
    ];
verifyTrue(testCase, all(startsWith(resolved, ...
    string(fullfile(testCase.TestData.appDir, 'functions')))));
verifyFalse(testCase, any(contains(resolved, 'legacy')));
end

function testRoiTraceExtractionMatchesDualMeanAndReportsLabels(testCase)
movie = reshape(uint16(1:(6 * 5 * 4)), 6, 5, 4);
mask = zeros(6, 5);
mask(2:3, 2:4) = 1;
mask(5:6, 1:2) = 3;
profile = roi_profile("calcium", [6 5]);
[actual, info] = extract_roi_traces(movie, mask, profile);

matrix = reshape(movie, 6 * 5, 4);
expected = [mean(matrix(mask(:) == 1, :), 1, 'double')', ...
    mean(matrix(mask(:) == 3, :), 1, 'double')'];
verifyEqual(testCase, actual, expected);
verifyEqual(testCase, info.labels, [1 3]);
verifyEqual(testCase, info.pixel_counts, [6 4]);
end

function testRoiTraceExtractionDoesNotPairOrAlterFrames(testCase)
movie = reshape(1:(4 * 4 * 7), 4, 4, 7);
mask = zeros(4, 4);
mask(1:2, 1:2) = 2;
profile = roi_profile("voltage", [4 4]);
[actual, info] = extract_roi_traces(movie, mask, profile);
verifySize(testCase, actual, [7 1]);
verifyEqual(testCase, info.labels, 2);
verifyEqual(testCase, info.pairing_rule, ...
    "none; channel pairing belongs to the section");
end

function testBackgroundMatchesFrozenFormula(testCase)
assumeEqual(testCase, exist('strel', 'file'), 2);
assumeEqual(testCase, exist('imdilate', 'file'), 2);
assumeEqual(testCase, exist('smooth', 'file'), 2);

ncols = 48;
nrows = 40;
nframes = 80;
[x, y] = ndgrid(1:ncols, 1:nrows);
movie = zeros(ncols, nrows, nframes, 'single');
for frame_idx = 1:nframes
    movie(:, :, frame_idx) = single(100 + 0.2 * x + 0.1 * y ...
        + sin(frame_idx / 9) + 3 * exp(-((x - 24).^2 + (y - 20).^2) / 20));
end
mask = zeros(ncols, nrows);
mask((x - 24).^2 + (y - 20).^2 <= 9) = 1;
profile = roi_profile("voltage", [ncols nrows]);

[actual, info] = remove_background(movie, mask, profile);
expected = frozen_background_formula(movie, mask, profile);
verifyEqual(testCase, actual.background_raw, expected.background_raw);
verifyEqual(testCase, actual.background_fitted, expected.background_fitted);
verifyEqual(testCase, actual.roi_minus_background_raw, ...
    expected.roi_minus_background_raw);
verifyEqual(testCase, actual.roi_minus_background_fitted, ...
    expected.roi_minus_background_fitted);
verifyEqual(testCase, actual.background_mask, expected.background_mask);
verifySize(testCase, actual.background_mask, [nrows ncols]);
verifyEqual(testCase, info.dual_bg_removed_stage, ...
    "roi_minus_background_fitted");
end

function testBackgroundRoleDoesNotChangeSharedFormula(testCase)
assumeEqual(testCase, exist('strel', 'file'), 2);
assumeEqual(testCase, exist('imdilate', 'file'), 2);
assumeEqual(testCase, exist('smooth', 'file'), 2);
movie = reshape(sin((1:(36 * 36 * 50)) / 31), 36, 36, 50);
[x, y] = ndgrid(1:36, 1:36);
mask = double((x - 18).^2 + (y - 18).^2 <= 9);
voltage = roi_profile("voltage", [36 36]);
voltage.roi.background.outer_distance = 14;
calcium = voltage;
calcium.role = "calcium";
actual_v = remove_background(movie, mask, voltage);
actual_c = remove_background(movie, mask, calcium);
verifyEqual(testCase, actual_v, actual_c);
end

function result = frozen_background_formula(movie, mask, profile)
[ncols, nrows, nframes] = size(movie);
matrix = reshape(movie, ncols * nrows, nframes);
labels = unique(mask(mask > 0), 'sorted');
background_raw = zeros(nframes, numel(labels));
background_fitted = zeros(nframes, numel(labels));
raw_corrected = zeros(nframes, numel(labels));
fitted_corrected = zeros(nframes, numel(labels));
background_mask = zeros(nrows, ncols);
mask_cells = cell(numel(labels), 1);
se1 = strel('disk', profile.roi.background.inner_distance);
se2 = strel('disk', profile.roi.background.outer_distance);
for idx = 1:numel(labels)
    roi = mask == labels(idx);
    annulus = imdilate(roi, se2) & ~imdilate(roi, se1);
    average_image = zeros(size(annulus));
    average_image(annulus) = mean(matrix(annulus(:), :), 2);
    sorted_pixels = sort(average_image(annulus));
    threshold_index = round(numel(sorted_pixels) ...
        * profile.roi.background.background_fraction);
    threshold = sorted_pixels(threshold_index);
    selected = (average_image <= threshold) & annulus;
    mask_cells{idx} = selected;
    background = mean(matrix(selected(:), :), 1);
    fitted = smooth(background, profile.roi.background.fit_span, ...
        char(profile.roi.background.fit_method))';
    roi_trace = mean(matrix(roi(:), :), 1);
    background_raw(:, idx) = background;
    background_fitted(:, idx) = fitted;
    raw_corrected(:, idx) = roi_trace - background;
    fitted_corrected(:, idx) = roi_trace - fitted;
end
for idx = 1:numel(labels)
    background_mask(mask_cells{idx}) = labels(idx);
end
result = struct( ...
    'background_raw', background_raw, ...
    'background_fitted', background_fitted, ...
    'roi_minus_background_raw', raw_corrected, ...
    'roi_minus_background_fitted', fitted_corrected, ...
    'background_mask', background_mask);
end

function profile = roi_profile(role, frame_size)
profile = struct( ...
    'role', role, ...
    'roi', struct( ...
        'frame_size', frame_size, ...
        'background', struct( ...
            'inner_distance', 2, ...
            'outer_distance', 32, ...
            'background_fraction', 0.05, ...
            'fit_span', 0.01, ...
            'fit_method', "rloess")));
end

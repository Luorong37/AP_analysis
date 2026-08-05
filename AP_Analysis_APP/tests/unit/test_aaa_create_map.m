function tests = test_aaa_create_map
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
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testVoltageMapMatchesFrozenDivisibleFormula(testCase)
movie = reshape(single(sin((1:(8*6*40))/17)), 8, 6, 40);
profile = map_profile("voltage", [8 6], 2);
actual = create_map(movie, profile);
expected = frozen_formula(movie, "voltage");
verifyEqual(testCase, actual, expected);
end

function testCalciumMapMatchesFrozenDivisibleFormula(testCase)
movie = reshape(single(cos((1:(8*6*40))/23)), 8, 6, 40);
profile = map_profile("calcium", [8 6], 2);
actual = create_map(movie, profile);
expected = frozen_formula(movie, "calcium");
verifyEqual(testCase, actual, expected);
end

function testNondivisibleGeometryFailsWithoutCropping(testCase)
movie = zeros(7, 6, 10);
profile = map_profile("voltage", [7 6], 4);
verifyError(testCase, @() create_map(movie, profile), ...
    'AAA:Functions:MapBinGeometryMismatch');
end

function profile = map_profile(role, frame_size, bin)
profile = struct('role',role,'roi',struct( ...
    'frame_size',frame_size,'map',struct('bin',bin)));
end

function map = frozen_formula(movie, role)
kernel = ones(3,3);
kernel(5) = 0;
for idx = 1:size(movie,3)
    movie(:,:,idx) = imfilter(movie(:,:,idx)./8, kernel, 'conv');
end
matrix = reshape(movie, size(movie,1)*size(movie,2), []);
map = zeros(size(matrix,1),1);
for idx = 1:size(matrix,1)
    if role == "voltage"
        corrected = detrend(double(matrix(idx,:)),1);
        point = find(abs(corrected) == max(abs(corrected)),1,'first');
        map(idx) = corrected(point);
    else
        corrected = detrend(double(matrix(idx,:)),2);
        map(idx) = std(corrected);
    end
end
map = reshape(map,size(movie,1),size(movie,2));
map(map>0) = map(map>0)-mean(map(map>0));
map(map<0) = map(map<0)-mean(map(map<0));
end

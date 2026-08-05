function tests = test_aaa_registration
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

function testManualPointsUseFrozenXYConvention(testCase)
movie = reshape(1:(6*5*8),6,5,8);
voltage = profile("voltage",[6 5]);
calcium = profile("calcium",[6 5]);
spec = struct('mode',"manual_points", ...
    'manual_points',[10 20; 7 15]);
[info, preview] = register_channels(movie,movie,voltage,calcium,spec);
verifyEqual(testCase, info.offset_xy, [3 5]);
verifyEqual(testCase, info.coordinate_rule, ...
    "voltage_position = calcium_position + offset_xy");
verifyEqual(testCase, size(preview.calcium_registered), [6 5]);
end

function testNoneUsesExplicitSavedOffset(testCase)
movie = zeros(4,7,3);
voltage = profile("voltage",[4 7]);
calcium = profile("calcium",[4 7]);
info = register_channels(movie,movie,voltage,calcium, ...
    struct('mode',"none",'reuse_offset',[-2 4]));
verifyEqual(testCase, info.offset_xy, [-2 4]);
verifyFalse(testCase, info.movie_transform ~= "none");
end

function testGeometryMismatchFails(testCase)
voltage_movie = zeros(4,7,3);
calcium_movie = zeros(7,4,3);
verifyError(testCase, @() register_channels( ...
    voltage_movie,calcium_movie,profile("voltage",[4 7]), ...
    profile("calcium",[7 4]),struct('mode',"none")), ...
    'AAA:Functions:RegistrationGeometryMismatch');
end

function value = profile(role, frame_size)
value = struct('role',role,'frame_rate',100, ...
    'transpose_before_analysis',false, ...
    'roi',struct('frame_size',frame_size));
end

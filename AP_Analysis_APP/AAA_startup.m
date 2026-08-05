function app_root = AAA_startup()
%AAA_STARTUP Add public scientific functions and the AAA backend.
%
% This intentionally avoids addpath(genpath(...)) so legacy snapshots and
% tests cannot shadow production functions/backend. The returned path is useful
% to App Designer callbacks and regression tests.

app_root = fileparts(mfilename('fullpath'));
backend_root = fullfile(app_root, 'backend');
functions_root = fullfile(app_root, 'functions');
normcorre_root = fullfile(app_root, '_external', 'NoRMCorre-0.1.1');

if ~contains(path, functions_root)
    addpath(functions_root, '-begin');
end
if ~contains(path, backend_root)
    addpath(backend_root, '-end');
end
if isfolder(normcorre_root) && ~contains(path, normcorre_root)
    addpath(normcorre_root, '-end');
end

fprintf('[AAA] Functions ready: %s\n', functions_root);
fprintf('[AAA] Backend ready: %s\n', backend_root);
end

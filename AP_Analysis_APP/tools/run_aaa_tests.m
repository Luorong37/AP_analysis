function [summary,report_file] = run_aaa_tests(test_targets)
%RUN_AAA_TESTS Run AAA tests and persist results before MATLAB process exit.
% This workstation's R2025a can fail later in DDUX shutdown. The report and
% AAA_TESTS_COMPLETE marker distinguish completed tests from that exit bug.
if nargin < 1 || isempty(test_targets)
    test_targets = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        'tests','unit');
end
app_root = fileparts(fileparts(mfilename('fullpath')));
addpath(app_root,'-begin');
AAA_startup();
suite = testsuite(test_targets,'IncludeSubfolders',true);
results = run(suite);
summary = table(string({results.Name})', [results.Passed]', ...
    [results.Failed]', [results.Incomplete]', seconds([results.Duration]'), ...
    'VariableNames',{'Name','Passed','Failed','Incomplete','DurationSeconds'});
report_dir = fullfile(app_root,'runtime','test_reports');
if ~isfolder(report_dir), mkdir(report_dir); end
tag = char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
report_file = string(fullfile(report_dir,['aaa_test_report_' tag '.mat']));
save(report_file,'summary');
disp(summary);
fprintf('[AAA_TESTS_COMPLETE] passed=%d failed=%d incomplete=%d report=%s\n', ...
    nnz(summary.Passed),nnz(summary.Failed),nnz(summary.Incomplete),report_file);
end

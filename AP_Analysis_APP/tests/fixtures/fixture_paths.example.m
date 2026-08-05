function paths = fixture_paths_example()
%FIXTURE_PATHS_EXAMPLE Copy to fixture_paths.m and fill local data paths.
%
% Fixed saved ROI and accepted-peak sources should be used for regression
% so manual interaction does not change the baseline.

paths = struct();
paths.ap_cycle = "";
paths.ap_source_results = "";
paths.dual_cycle = "";
paths.dual_source_results = "";
paths.dual_record = "";
end

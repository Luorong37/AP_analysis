function receipt = AP_analysis3(request,runtime)
%AP_ANALYSIS3 Internal compatibility name for the shared single runtime.
if nargin < 2, runtime = struct(); end
receipt = aaa.workflows.run_single(request,runtime);
end

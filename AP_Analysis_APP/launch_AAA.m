function app = launch_AAA()
%LAUNCH_AAA Start the AP Analysis APP from any current MATLAB folder.

project_root = fileparts(mfilename('fullpath'));
if ~contains(path, project_root)
    addpath(project_root, '-begin');
end
AAA_startup();
app = AP_Analysis_APP;
if nargout == 0
    clear app
end
end

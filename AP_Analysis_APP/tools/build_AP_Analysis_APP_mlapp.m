function target = build_AP_Analysis_APP_mlapp()
%BUILD_AP_ANALYSIS_APP_MLAPP Build the editable App Designer container.
%
% The class source is kept as text so this builder can be reviewed and
% regenerated.  Scientific code is never embedded in the App.

tool_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(tool_dir);
source_file = fullfile(tool_dir, 'source', 'AP_Analysis_APP.m');
target = fullfile(project_root, 'AP_Analysis_APP.mlapp');

code = fileread(source_file);

fig = uifigure('Visible', 'off', 'Name', 'AP Analysis APP (AAA)');
fig.Position = [40 40 1500 900];
grid = uigridlayout(fig, [2 1]);
grid.RowHeight = {42, '1x'};
title_label = uilabel(grid, 'Text', 'AP Analysis APP (AAA)');
title_label.FontSize = 18;
title_label.FontWeight = 'bold';
note_label = uilabel(grid, 'Text', ...
    ['The runtime interface is generated from the AAA request, section, ' ...
     'and parameter catalogs.']);
note_label.HorizontalAlignment = 'center';

serializer = appdesigner.internal.serialization.MLAPPSerializer(target, fig);
serializer.OverwriteTargetFile = true;
serializer.MatlabCodeText = code;
serializer.ClassName = 'AP_Analysis_APP';
serializer.Metadata = appdesigner.internal.model.MetadataModel();
serializer.save();

delete(note_label);
delete(title_label);
delete(grid);
delete(fig);

fprintf('[AAA] App Designer file written: %s\n', target);
end

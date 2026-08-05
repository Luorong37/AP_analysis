function fig = plot_trace_overview(channel_results,profiles)
%PLOT_TRACE_OVERVIEW Render a mode-neutral processed-trace overview.
if ~isstruct(channel_results) || ~isstruct(profiles) ...
        || numel(channel_results) ~= numel(profiles)
    error('AAA:Functions:InvalidTraceOverviewInput', ...
        'Results and profiles must be paired struct arrays.');
end
fig = figure('Color','w','Name','AAA Trace Overview');
layout = tiledlayout(fig,numel(profiles),1,'TileSpacing','compact');
for idx = 1:numel(profiles)
    ax = nexttile(layout);
    [data,stage] = resolve_trace_stage(channel_results(idx), ...
        {'sensitivity_smoothed','sensitivity','bleach_removed','raw'});
    time = (0:size(data,1)-1)' / double(profiles(idx).frame_rate);
    plot(ax,time,data,'LineWidth',.8);
    grid(ax,'on'); xlabel(ax,'Time (s)');
    ylabel(ax,string(profiles(idx).role));
    title(ax,sprintf('%s — %s',profiles(idx).role,stage), ...
        'Interpreter','none');
end
end

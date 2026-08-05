function [peaks_index, peaks_amplitude, peaks_polarity] = ...
        peak_table_to_cells(peak_table, nrois, fallback_polarity)
%PEAK_TABLE_TO_CELLS Convert accepted rows to the frozen Dual cell layout.
% Inputs:
%   peak_table Table containing roi, index, amplitude_sensitivity, polarity,
%              and status variables. Indices are one-based frames.
%   nrois Nonnegative integer output ROI count.
%   fallback_polarity Optional 1-by-ROI cell array used for empty ROIs.
% Outputs:
%   Three 1-by-ROI cell arrays ordered by ROI and ascending frame index.

if nargin < 3
    fallback_polarity = cell(1, nrois);
end
if ~istable(peak_table)
    error('AAA:Algorithms:InvalidPeakTable', 'peak_table must be a table.');
end
required = {'roi','index','amplitude_sensitivity','polarity','status'};
if ~all(ismember(required, peak_table.Properties.VariableNames))
    error('AAA:Algorithms:InvalidPeakTable', ...
        'peak_table is missing required variables.');
end
validateattributes(nrois, {'numeric'}, ...
    {'scalar','integer','nonnegative','finite'});

peaks_index = cell(1, nrois);
peaks_amplitude = cell(1, nrois);
peaks_polarity = cell(1, nrois);
for roi_idx = 1:nrois
    rows = peak_table(peak_table.roi == roi_idx ...
        & peak_table.status == "accepted", :);
    if height(rows) > 0
        [~, order] = sort(rows.index);
        rows = rows(order, :);
        peaks_index{roi_idx} = rows.index;
        peaks_amplitude{roi_idx} = rows.amplitude_sensitivity;
        peaks_polarity{roi_idx} = rows.polarity(1);
    else
        peaks_index{roi_idx} = [];
        peaks_amplitude{roi_idx} = [];
        if numel(fallback_polarity) >= roi_idx ...
                && ~isempty(fallback_polarity{roi_idx})
            peaks_polarity{roi_idx} = fallback_polarity{roi_idx};
        else
            peaks_polarity{roi_idx} = 1;
        end
    end
end
end

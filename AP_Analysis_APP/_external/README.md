# External runtime dependencies

`NoRMCorre-0.1.1` is copied verbatim from the installed analysis runtime so
AAA does not depend on a previously configured MATLAB search path.  Its
upstream license is retained in the dependency directory.

Project-owned compatibility helpers (`calculate_FWHM.m`, `offset_plot.m`,
and `array2tif.m`) are frozen beside the legacy workflows under
`legacy/snapshot_2026-07-30` because those workflows call them directly.

# AAA migration status

## Production path now migrated

- `AP_Analysis_APP.mlapp` is the public entry and calls
  `aaa.run_analysis3(request,runtime)`.
- AP/Dual Single share `aaa.workflows.run_single`, `create_context`, the
  ordered section registry, and the same `run/reuse/skip` runtime.
- Production no longer executes the frozen scripts or a base-workspace
  compatibility bridge.
- Scientific/data functions are ordinary external `.m` files under
  `functions`; package code under `backend/+aaa` is limited to schema,
  orchestration, handlers, logging, and persistence.
- Input, rigid motion, Dual registration, ROI, trace processing, peaks,
  AP/Dual automatic visual-stim analysis, the frozen Dual overview/overlap/
  tuning/frequency figure bundles, Dual FFT/CWT, and export have handlers.
- Every completed/reused section now writes one `AAA.result.v1` container
  under `results/`; trace-stage provenance, including `info.parent_results`,
  is preserved while historical main-container filenames are no longer written.
- `load_analysis_results.m` is the single old-format compatibility boundary.
  Reuse discovery, section reuse, and Record loading consume its normalized
  `results.voltage`/`results.calcium` schema and conversion report.
- ROI reuse does not recompute sensitivity maps. Record runners delegate
  per-Cycle work to the same Single workflow.
- AP/Dual trace reuse now auto-resolves blank Cycle/Record sources, restores
  result containers without loading movies, and records each chosen source
  in App/session/Record logs and job plans.
- The App header exposes a reuse source/search-root field plus Browse and
  Detect controls. Record-shaped input folders automatically select Record
  scope; nested `Rec*/Cycle*` folders are included.
- Run Log is fixed below the main tabs. Input/path/control edits are logged,
  and directory preflight inventories every Cycle and saved result before Run.

## Explicitly not migrated in this phase

- AP nonrigid motion and VolPy are temporarily removed from the active
  request/UI/runtime.
- AP event windows/FWHM and AP statistics. `run_single` records these in
  `excluded_sections` and forces them to `skip` even if an old preset says
  `run`.
- Real-data numerical certification of manifest/log variants and every
  restored presentation figure remains pending fixture selection.

## Verification completed

- Clean-path location tests confirm scientific functions resolve from
  `AP_Analysis_APP/functions`, not a package or legacy folder.
- Unit tests cover input/motion, registration, map parity, ROI atoms and
  handler, trace atoms and handler, peaks, comparison, stimulus discovery and
  metrics, tuning, output renderers, frequency analysis, section runtime,
  schema, IO, and manifest creation.
- Dedicated reuse tests cover newest-completed selection, exact per-Cycle
  Record matching, no-movie result restoration, and creation of a new AP
  reuse output folder without any movie file.
- A real historical AP result folder with nonstandard Single filenames was
  loaded by MAT contents at 400 Hz, then trace-reuse converted into standalone
  `results/roi.mat` plus `results/trace.mat` without loading a movie.
- Dual registry: all 11 sections have handlers.
- AP registry: all sections except AP events/statistics have handlers,
  including the shared visual-stim handler.

## Approval boundary

Numerical certification still needs user-selected AP Single, Dual Single,
and Dual Record fixture directories. Structural and synthetic tests do not
certify exact equality for every real acquisition/output figure.

MATLAB R2025a on this workstation also intermittently crashes after all tests
have passed, in `mwddux_matlab.dll` during process exit. Test tables are
captured before exit; this is tracked separately from AAA failures.

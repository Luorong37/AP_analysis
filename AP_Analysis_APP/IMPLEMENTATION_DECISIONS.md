# AAA implementation decisions

This file records implementation choices made while constructing the first
review version.  It is deliberately explicit so that UI convenience cannot
silently change algorithm meaning.

1. The existing parent-folder analysis files are not edited.  The migration
   baseline is copied into `legacy/snapshot_2026-07-30` with source hashes.
2. Scientific/data functions are ordinary external `.m` files under
   `functions`; only control infrastructure is namespaced under `aaa`.
3. AAA is the public entry. Public same-name AP/Dual/Record facades and
   base-workspace scientific execution were removed.
4. Single and Record share the same request schema.  Record is orchestration,
   not a second copy of Single algorithms.
5. A reference cycle is represented as a planned job.  The planner must not
   schedule the same full reference job twice when it is also included among
   ordinary cycles.
6. `ap_statistics` initially groups AP trend/average,
   sequence-density, and ISI/phenotype analyses.  It can later be split by
   adding catalog entries without changing the dispatcher contract.
7. AP Record v1 dispatches Single cycles only.  No new AP Record averaging
   algorithm is inferred from the ignored legacy batch file.
8. Dual Record may run `Dual_rec_analysis3` once after all successful cycle
   jobs, using exact result paths from the job manifest.
9. Section actions are tri-state (`run`, `reuse`, `skip`).  Presets only
   initialize this table; an edited table becomes the explicit custom plan.
10. Every logical default/policy is exposed as a checkbox.  Nonlogical
    defaults are exposed in Basic, Advanced, or Expert tables.  Parameters
    may be conditionally relevant, but remain discoverable.
11. The App does not show a fake percentage or offer unsafe mid-algorithm
    cancellation.  It shows section/cycle events emitted at real boundaries.
12. The App's return receipt is lightweight.  Scientific results remain in
    the target folder and preserve the legacy file/data structure.
13. `analysis.log` and `analysis_manifest.mat` are additive traceability
    outputs and do not replace legacy result files.
14. Empty output paths, run names, and optional reuse paths retain the legacy
    naming/fallback rule; the selected fallback is recorded in log/manifest.
15. Regression parity cannot be certified without user-provided AP Single,
    Dual Single, and Dual Record fixtures.  Until then, structural smoke tests
    and static checks are evidence of wiring only, not numerical equality.
16. The `.mlapp` UI is catalog-driven: logical rows use native checkbox
    columns, section actions use categorical dropdowns, and nonlogical values
    use editable Basic/Advanced/Expert tables.  This avoids one static widget
    per parameter and makes additions/removals schema-driven.
17. The App-generated request is the only public control source. Backend
    compatibility workflow names remain only so Record scheduling has a
    stable internal call target.
18. Public values that the current backend cannot apply are visible with
    their real defaults and migration rule, but read-only. Validation rejects
    an attempted override before scientific execution. Parameters already
    consumed by migrated atoms (rigid motion, Dual noise, and frequency) are
    editable.
19. AP and Dual minimum peak-distance defaults are stored as empty overrides,
    preserving their frame-rate-dependent legacy formulas.
20. `skip_existing` trusts a completed AAA manifest when one exists; failed
    or pending tracked runs are not reused.  Old result folders without a
    manifest remain compatible and may still be selected by legacy markers.
21. NoRMCorre 0.1.1 and directly called helper files are vendored inside the
    isolated project.  `AAA_startup` adds explicit roots rather than relying
    on the parent repository's saved MATLAB path.
22. `analysis.log` is mandatory and remains a public read-only checkbox set
    to true.  `analysis_manifest.mat` remains default-on but optional, matching
    the approved optional-manifest design.
23. VolPy controls and execution are temporarily removed from the active AAA
    request/UI/runtime. Frozen source remains only as an audit snapshot.
24. ROI reuse skips activity-map computation because maps guide drawing but
    do not affect extraction from an already accepted mask. This removes a
    repeated Record-Cycle computation; the manifest/action still records ROI
    reuse explicitly.
25. Registration offsets are consistently `[x y]`; earlier UI text saying
    `[row column]` was documentation error and was corrected.
26. AP events/statistics remain excluded rather than approximated. AP
    nonrigid and its inactive public controls are temporarily removed; rigid
    NoRMCorre remains active.
27. AP and Dual share manifest/log/sync stimulus discovery and the frozen
    baseline/stim mean, peak, AUC, and paired-test formulas. Voltage response
    is public: `auto` uses accepted peak count when peak analysis is present
    and otherwise falls back to mean sensitivity, with the effective method
    recorded in results and logs.
28. The frozen map centering order is preserved literally: positive pixels
    are centered first, then the negative mask is recomputed from the changed
    map and centered. Caching both masks beforehand changes the numbers.
29. MATLAB R2025a process-exit crashes in `mwddux_matlab.dll` are treated as
    environment failures only after the test framework has printed passing
    results. Startup tests use at least a 120-second external timeout because
    cold startup on this workstation takes roughly 45-58 seconds.
30. `reuse_trace` and `analysis_only` restore saved result containers and a
    derived time axis without loading movie pixels. `reuse_motion` still
    loads the movie because applying saved shifts is a pixel operation.
31. A blank reuse path triggers automatic detection by default. An explicit
    path may be either one result folder or a Cycle/Record search root. A
    direct compatible result folder wins; otherwise the newest compatible
    completed candidate is selected by manifest `completed_at`, falling back
    to file time for legacy results.
32. Record reuse is resolved independently for each exact Cycle path
    component. Direct and nested `Rec*/Cycle*` layouts are supported. Search
    is narrowed to the target Cycle before result anchors are scanned so the
    same Record tree is not recursively traversed once per candidate file.
33. A tracked AAA result whose manifest is not `completed` is never reused.
    A legacy result without a manifest remains eligible when its required
    section files are complete.
34. No trace values are interpolated, padded, truncated, renormalized, or
    recalculated during trace reuse. Frame count comes from saved
    `trace_results.raw.data`. The saved trace frame rate is authoritative,
    replaces the UI/request value for downstream time conversion, and both
    values are logged when they differ.
35. The new reuse output copies only frozen upstream MAT artifacts needed for
    provenance and downstream readers (ROI/map/registration/background/
    bleach/metric/dual container where present). It does not copy movies or
    silently recompute missing upstream products.
36. Choosing a folder named `Cycle*` selects Single scope; choosing a folder
    containing direct `Cycle*` or nested `Rec*/Cycle*` folders selects Record
    scope. Automatic scope changes preserve the selected workflow and
    mode-specific parameters, initialize only the new scope controls, retain
    the input path, and log the old/new scope.
37. Input owns movie loading only and permits `run` or `skip`. `reuse_trace`
    and `analysis_only` use `input=skip, trace=reuse`; Trace alone restores
    result containers, metadata, time axes, and upstream provenance MAT files.
38. Input preflight reads directory entries and MAT variable directories but
    never movie arrays. A movie-dependent plan is blocked before execution if
    any selected Cycle is empty or lacks a standard loadable source. A
    trace-reuse plan may proceed without movie input only when every Cycle has
    a compatible completed result.
39. A nonempty manual Reuse source is authoritative. Automatic candidates are
    still scanned and logged for comparison but never replace it. An invalid
    manual path blocks execution rather than falling back automatically.
40. The App logs old/new values for paths, mode, scope, preset, sections, and
    public parameters. Its Run Log remains in the main bottom panel; the
    parameter and section tabs remain in the center panel.
41. Result filenames are discovery hints, not compatibility evidence.
    `analysis_manifest.mat` output references and `*_results.mat` locate
    candidate folders; every selected trace must then expose the shared
    `trace_results` stages with consistent matrix sizes, an unambiguous role,
    and a positive saved frame rate. Legacy folders without a manifest are
    labeled semantic legacy results and remain read-only.
42. Movie discovery prefers Cycle manifest metadata, resolves a missing
    logical TIFF into its `NAME_stackNN` parts, then falls back to matching
    camera folders. Preflight and the input section call the same resolver.
    They do not open movie pixels during discovery.
43. Manifest schema 1.1 records, per channel, the logical movie path,
    resolved source, ordered physical files, split flag and part count,
    resolution method, frame rate/count, and reuse result provenance. A
    historical source that has moved remains reusable; unavailable original
    movie-part provenance is left empty rather than guessed.
44. Reuse discovery is bounded by public Advanced controls (default depth 4
    and 500 MAT anchors per Cycle). Hitting a limit disables automatic
    selection and requires a manual source instead of silently choosing from
    an incomplete scan.
45. AP Record `cycle_receipts` use the complete single-run orchestration
    receipt field set and field order for completed, failed, and skipped
    cycles. This is a control-schema correction only; scientific result
    files and numerical behavior are unchanged.
46. The main Run Log starts at 320 pixels high and shares a draggable
    horizontal divider with the Sections/parameter tabs. Dragging is clamped
    to a 160-pixel minimum Log and a 220-pixel minimum tabs area. The chosen
    size lasts for the current App session only and is not stored in request
    or analysis manifests.
47. Record-level analysis remains outside the per-Cycle section table. The
    existing `record.record_average` switch is authoritative; the dependent
    `record_average_only` switch selects existing Cycle results and is invalid
    unless Record analysis is enabled. Both switches live in the Request panel
    and are hidden from the duplicate Basic/Advanced parameter tables.
48. `Dual_rec_analysis3.m` is frozen and was not edited. The independent
    public `Record_analysis3.m` implements the role-aware AP/Dual Record
    interface and calls external functional helpers. AP Record scheduling now
    calls this entry; the production Dual scheduler continues calling the
    frozen Dual implementation so its established files and figures do not
    regress before numerical-equivalence approval.
49. Record trace averaging inherits the Dual rules: flash aligns to first
    flash onset, grating aligns to recording start, only the common overlapping
    relative-time interval is retained, interpolation is linear, and saved
    frame-rate or ROI-count mismatches are errors. Saved result frame rates are
    authoritative and are logged per Cycle/role.
50. Record grating curves use the union of directions modulo 360; missing
    directions remain NaN. Both the mean of per-Cycle metrics and metrics
    recomputed from the Cycle-mean curve are saved. Peak-count and
    mean-sensitivity voltage responses are never pooled; each records its own
    included Cycle list. A legacy voltage tuning without method metadata is
    classified as peak count to preserve the established Dual interpretation.
51. Result filenames are hints during Record loading. A direct result folder
    is inspected for scalar channel/stim payload variables, permitting
    nonstandard MAT filenames. Record manifests save included result folders,
    saved frame rates, source movie-part paths, and part counts. AP outputs use
    `AP_analysis3_rec_average`; `record_average_only` never loads movie pixels.
52. AP flash Record FFT/CWT is explicitly controlled by public AP Record
    parameters. Its inherited defaults are enabled, 0.5-100 Hz, `amor`, and
    12 voices per octave; the maximum is clamped below the saved Nyquist
    frequency. Grating Records do not run the flash frequency module.
53. Record source-resolution provenance is stored as one cell per Cycle because
    resolver reports may expose different diagnostic fields. This changes only
    traceability storage; source selection and numerical analysis are unchanged.
54. Every App push button writes and renders an immediate `[UI]` click entry.
    Dialog cancellation, blocked actions, errors, and successful output opening
    are also logged so a modal dialog or long preflight cannot look unresponsive.
55. Input/reuse preflight results are cached for the current App session using
    the complete serialized request as the cache key. Any request change causes
    a new inspection; the explicit Detect Reuse source button always refreshes.
    External filesystem changes require that explicit refresh or a request edit.
56. Preset, Section Action, Record switches, paths, and public parameter edits
    now update request state and Log immediately, then invalidate the session
    preflight record without scanning result folders. The log records old/new
    values where available and the control that invalidated validation.
57. An Input path edit performs only a quick structural check: existence,
    inferred Single/Record scope, Cycle count/emptiness, and standard movie
    source evidence. It does not inspect saved result MAT files or decide reuse
    compatibility; those decisions are deferred to Detect, Validate, or Run.
58. Detect and Validate create an authoritative session preflight record. Run
    reuses a matching passed or failed record; a failed record blocks without
    rescanning. Without a matching record, Run performs one authoritative
    preflight. Selected source paths receive only a cheap existence check before
    reuse; a missing path forces a new preflight.
59. The authoritative record keeps the user request unchanged and separately
    stores the resolved execution request, including automatic Single-Cycle
    reuse source and saved frame-rate authority. It is passed through runtime so
    dispatcher, Single, AP Record, and Dual Record backends skip matching deep
    scans. Backend calls without this record retain their standalone preflight.
    Record child matching uses the exact normalized Cycle path and never crosses
    Cycle boundaries.
60. Runtime progress events are emitted at functional phase boundaries and per
    role, not per frame. ROI selection and manual peak gating explicitly report
    `waiting`; Record Cycle-result loading reports `n/N`. These events change
    only logging/UI status and do not change scientific calculations or files.
61. Record loading accepts both the historical complete channel-result schema
    (`movie_info` plus `trace_results`) and the shared AP trace-stage schema.
    For the latter, `raw.frame_rate` is the authoritative saved rate and the
    stage data row count is the frame count. A manifest rate, when present and
    valid, is cross-checked at 1e-9 relative/absolute tolerance; a mismatch is
    an error. Manifest source paths and split-file provenance are copied into
    the normalized in-memory `movie_info`. UI/request frame rate is never used
    as a fallback, and saved trace values are not recalculated.
62. The App header lamp is green whenever AAA is idle and red only while the
    analysis dispatcher is executing. Validate and Detect Reuse remain green.
    Existing Run cleanup returns the lamp to green after both success and error;
    no persistent error color or third state is introduced.
63. Random-grating Record visualization now reads saved per-Cycle sensitivity
    stages and saved baseline/stimulus frame windows. Trials include the full
    baseline through stimulus end, align stimulus onset to zero, are never
    interpolated, and use relative-frame NaN padding only in heatmaps. Directions
    are sorted after modulo 360; within one direction, Trials retain Cycle order
    and original Trial order.
64. Trial polarity first uses the per-role polarity saved with stimulus
    sensitivity metrics. Historical results without it fall back to the explicit
    Record control polarity; no polarity is inferred during Record rendering.
    Stack spacing inherits the Dual display rule of six times median Trial
    standard deviation, with median range (minimum one) only as a display fallback.
    Sorted directions use the same fixed HSV color order in Trial stacks and
    continuous Cycle traces; color has no analytical meaning.
65. Detailed Trial figures are generated only for saved random-grating records,
    identified by `windows.stim_type` or the historical
    `random_grating_tuning` marker. Ordinary grating Records keep tuning outputs
    without the large Trial figure set. Heatmaps use one symmetric color limit
    per role from all displayed Trial sensitivity values. Per-Cycle DSI/OSI
    figures consume saved `metrics_per_cycle` and do not recompute metrics.
66. Random-grating Record outputs use external role-aware rendering code and
    save FIG plus 150-dpi PNG bundles under `9_rg_stack`, `10_rg_hm`,
    `11_rg_trace`, and `12_rg_cycle_metrics`. AP emits voltage only; Dual emits
    its available roles. Calcium peak-raster algorithms and the frozen Dual
    scheduler remain unchanged.
67. A multi-role detailed Record render requires equal ROI counts across roles
    because one combined ROI/Cycle figure uses one ROI index for every role.
    This is validated before any detailed figure is written. AP has one role
    and is unaffected; the rule inherits the paired-ROI assumption of the
    historical Dual Record plots rather than padding or dropping ROIs.
68. Structured Cycle output now uses only `results/<section>.mat`; every file
    contains one `result` variable with schema `AAA.result.v1`. The public
    loaded container exposes `results.voltage`, optional `results.calcium`,
    `results.stim`, and `results.record` directly, without `channel`,
    `role_results`, or `artifact` container levels.
69. `save_analysis_results(ctx,section,outcome)` is the sole section-result
    writer and is called by `run_plan` after every successful or reused
    handler. Skip writes nothing. PNG/TIFF/FIG/CSV/log files remain separate
    `output_files`; motion shift, parameter, and QC MAT files remain separate
    because motion reuse/application consumes them as algorithm files.
70. Raw ROI traces are stored in `results/roi.mat` and removed from
    `results/trace.mat`. When trace-only reuse skips ROI execution, the saver
    creates one reused `roi.mat` dependency from the loaded raw traces and
    masks before writing `trace.mat`; raw is therefore not duplicated and the
    destination remains independently readable.
71. `load_analysis_results.m` is the only historical AP/Dual result adapter.
    Filenames are discovery hints only. If several semantic containers map to
    one role, the loader selects the container with more complete trace/peak/
    movie metadata and uses newest file time only as a tie-break; the choice
    and rejected ambiguity are written to `report.warnings`.
72. Every compatibility report records original and absolute source paths,
    absolute and relative source files, detected schema, adapter, mapped and
    missing fields, warnings, load time, and saved frame-rate source. Reuse
    logs this report, saves it in result provenance, and adds one manifest
    history entry per source/adapter pair.
73. Saved movie/trace frame rate remains authoritative. A legacy manifest
    frame rate is copied only when trace metadata lacks one; disagreement is
    an error at the existing 1e-9 absolute/relative tolerance. Source movie
    parts and split-file metadata are copied from the manifest without reading
    movie pixels.
74. Record averages now write `results/record.mat` through the same saver.
    Renderer-only MAT duplicates (`rec_avg.mat`, per-stage MAT, grating MAT,
    and Record frequency MAT) are removed; their scientific values live in
    `result.data`, while FIG/PNG/CSV remain `output_files`. Both AP and Dual
    Record schedulers call the shared external `Record_analysis3.m`; the old
    `Dual_rec_analysis3.m` remains untouched as a reference.
75. Dual overview rendering now consumes in-memory/shared maps, ROI masks,
    trace stages, and background data instead of reloading historical
    intermediate MAT files. The derived heatmap MAT duplicate is removed;
    figure calculations, polarities, scaling, and plotted values are unchanged.
76. No event alignment, averaging, baseline, response-window, smoothing,
    filtering, interpolation, polarity, peak, or tuning algorithm was changed
    by the unified result migration.
77. Dual Record `skip_existing` and `record_average_only` now recognize a
    Cycle result from a completed `analysis_manifest.mat` containing a
    completed/reused `trace` result record. A completed motion-only manifest
    is therefore insufficient. Historical `-1_explicit_dual_results.mat`
    remains a read-only fallback, and all compatible candidates share the
    same newest-timestamp selection. The `analysis_only` search additionally
    requires the manifest request preset (or legacy `dual_info.mat` preset)
    to be `analysis_only`.
78. Section handlers now expose the direct order `Load inputs -> Run algorithm
    -> Return results`; no generic wrapper or renamed atom was introduced.
    Existing external scientific function names are preserved and their
    headers document actual array dimensions, units, required bound-profile or
    spec fields, and outputs. Channel profiles retain only resolved role,
    frame-rate/orientation, motion, ROI, trace, and peak settings; the former
    `profile.parameters.common/mode` copy of the request was removed. The Dual
    raw-pair source rule remains in the input handler and reads its public
    source-discovery settings directly from `ctx.request.params.dual`. This
    boundary change does not alter any scientific formula, numeric default,
    reuse decision, alignment rule, or saved-result schema.
79. AAA release identity is stored in the root `VERSION` file and documented
    in `CHANGELOG.md`. New manifests record that App version, repository HEAD
    commit/branch when directly readable, MATLAB release, NoRMCorre version,
    and an ordered SHA-256 digest of production source. Runtime provenance does
    not invoke Git because external Git processes have caused blocking on the
    target workstation. Worktree cleanliness is therefore honestly recorded
    as `not_checked`; the digest, rather than an inferred clean flag, identifies
    the actual local source. Source snapshot copying is publicly visible but
    disabled because it was not implemented and would duplicate the manifest
    provenance without a single authoritative policy. No numerical behavior,
    file naming, result schema, or reuse decision changes with this policy.

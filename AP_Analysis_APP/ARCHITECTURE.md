# AAA architecture

```text
AP_Analysis_APP.mlapp
  -> request + explicit section plan
  -> aaa.run_analysis3
       -> Single: aaa.workflows.run_single
       -> Record: cycle planner -> run_single once per planned Cycle
            -> aaa.sections.run_plan
                 -> thin handler
                      -> Load inputs (request/reuse/files)
                      -> Run algorithm (ordinary functions/*.m)
                      -> Return results (ctx/outcome/log)
                 -> save_analysis_results(ctx, section, outcome)
```

The `.mlapp` owns UI state only. It does not contain movie, ROI, trace,
peak, comparison, or frequency algorithms. The lightweight receipt contains
status and paths; scientific arrays stay in result files.

## Boundaries

- `functions/*.m`: external, directly testable functional scientific/data
  units. Existing function names are retained. Each function header records
  its real dimensions, units, required profile/spec fields, and output format.
  A bound channel profile is passed when role-dependent behavior is required,
  so role and parameters cannot drift between duplicated arguments. A profile
  contains effective algorithm settings but never retains the complete request.
- `backend/+aaa/+sections`: thin action handlers; all structured persistence
  crosses the external `save_analysis_results.m` boundary. Handler main bodies
  use the visible order `Load inputs -> Run algorithm -> Return results`;
  reuse loads and returns without calling a calculation atom.
- `backend/+aaa/+schema`: the only request/default/section catalog.
- `backend/+aaa/+io`: events, analysis log, and manifest.
- `legacy/snapshot_2026-07-30`: read-only behavioral reference, never a
  production dependency.

Manifest creation calls `aaa.io.source_provenance`. It reads repository HEAD
metadata without starting Git and hashes the ordered production source bundle
with SHA-256. Consequently a clean release is identified by its version/tag
and commit, while a run made from local edits is still identified by the exact
source digest. Tests, runtime output, and the frozen legacy snapshot are not
part of that digest; the vendored dependency is represented by its pinned
version and `DEPENDENCY_HASHES.txt`.

AP and Dual use the same input, motion, ROI, trace, peak, stimulus,
visualization, and export handlers. Single and Record never contain duplicate
scientific paths. Dual-only registration/comparison/frequency handlers call
ordinary external functions.

## Data and coordinate rules

- `profile.role` is the sole role authority.
- Trace stages retain explicit ordered `parent_results`.
- Dual ROI offset is `[x y]` with
  `voltage_position = calcium_position + offset_xy`.
- Movies are not moved by registration.
- Different channel frame counts are not truncated, padded, or interpolated.
- Map geometry must be exactly divisible by the configured map bin; there is
  no hidden crop/pad fallback.
- Strict reuse loads complete section results and fails when they are
  incomplete; it does not silently recompute from an earlier stage.
- `resolve_reuse_source` is the shared Single/Record AP/Dual locator. It
  validates the active section plan and completed manifest, matches Record
  candidates by exact Cycle component, and selects the newest valid result.
- Trace/analysis reuse has no `movie_3d` in its context. Only saved result
  containers, metadata, and time axes are restored. Motion reuse is the
  deliberate exception because it must apply shifts to frames.
- Input has no reuse meaning: it is `run` only for movie ingestion and `skip`
  for trace/analysis reuse. Directory preflight blocks invalid movie inputs
  and unresolved per-Cycle reuse before a workflow creates partial outputs.

## Result layout

```text
analysis_manifest.mat
results/
  input.mat
  motion.mat
  registration.mat
  roi.mat
  trace.mat
  peak.mat
  stim.mat
  comparison.mat
  frequency.mat
  visualization.mat
  export.mat
  record.mat          # Record-average workflow only
```

Every file stores one variable named `result` using `AAA.result.v1`.
`load_analysis_results.m` is the only historical AP/Dual compatibility
boundary and returns role fields directly as `results.voltage` and, when
present, `results.calcium`. PNG/TIFF/FIG/CSV/log files and motion shift/QC
files are referenced through `output_files`; they are not result containers.

# AP_Analysis_APP (AAA)

AAA is the isolated App Designer project for AP and dual-camera analysis.
The original files in the parent `AP_analysis` directory are not modified.

## Design rules

1. **Extensible:** AP/Dual and Single/Record share one request, parameter
   catalog, section catalog, validator, event protocol, and manifest format.
2. **Public:** every preset/default decision is represented in
   `aaa.schema.parameter_catalog` and is available to the App at Basic,
   Advanced, or Expert level.
3. **Simple:** scientific algorithms are shared only when their numerical
   meaning is verified to be the same. Single-cycle workflows own all
   per-cycle computation; record workflows only plan and dispatch jobs.
4. **Traceable:** each run writes `analysis.log`, `analysis_manifest.mat`,
   and one `results/<section>.mat` file per completed/reused section.
5. **Compatible:** new AAA writes only `AAA.result.v1` structured results;
   `load_analysis_results.m` reads both this schema and historical AP/Dual
   containers and records the adapter, absolute source files, mappings,
   warnings, and saved frame-rate authority.
6. **Explicit function contracts:** section handlers visibly load inputs, run
   existing-name algorithms, and return results. Scientific functions receive
   arrays plus a narrow bound profile/spec, never the App request or context.

The isolated project keeps the frozen source snapshot only as an auditable
numerical reference and vendors NoRMCorre 0.1.1 under `_external`.
Production execution does not run the frozen scripts. `AAA_startup` adds the
ordinary `functions` folder, the namespaced backend, and this explicit
dependency root; it does not scan the parent repository with `genpath`.

## Start

```matlab
app = launch_AAA;
```

Alternatively, double-click `AP_Analysis_APP.mlapp` in MATLAB and press
Run.  The `.mlapp` adds only the namespaced backend through `AAA_startup`.

The reviewable App class source is `tools/source/AP_Analysis_APP.m`. Rebuild
the `.mlapp` with `tools/build_AP_Analysis_APP_mlapp.m` after changing that
source; parameter/section additions normally require only catalog changes.

Backend-only calls use the package namespace:

```matlab
request = aaa.schema.default_request("dual", "single");
receipt = aaa.run_analysis3(request, struct());
```

The receipt contains status and paths only. Scientific results remain in
the target result directory.

For tests, temporarily add `tools` and call `run_aaa_tests`. It writes a MAT
summary and prints `[AAA_TESTS_COMPLETE]` before process exit, so a later
MATLAB DDUX shutdown fault cannot be mistaken for a failed test suite.

The App is the public entry. It constructs an explicit request and section
plan, then calls `aaa.run_analysis3`. Compatibility workflow names remain
inside the backend only for Record scheduling; there are no public
same-name script/function facades in the project root.

## Version management

The human-readable App version is stored in `VERSION`; release notes are kept
in `CHANGELOG.md`. Every new manifest records that version, the repository
commit readable from `.git/HEAD`, the MATLAB release, the pinned NoRMCorre
version, and a deterministic SHA-256 digest of the production source bundle.
The App does not launch a Git process during analysis, so provenance capture
cannot block a run on a slow Git command. `git_worktree_status` is therefore
explicitly `not_checked`; the source digest identifies the actual files used,
including uncommitted edits.

Suggested release rhythm:

1. Make and test one coherent change on a short-lived branch.
2. Update `CHANGELOG.md`; increment `VERSION` when publishing a release.
3. Commit source and documentation, but never `runtime/` outputs.
4. Tag tested releases as `aaa-vX.Y.Z` and push the branch and tag to GitHub.

Development before AP event/FWHM/statistics parity remains in the `0.x`
series. Use a patch increment for fixes, a minor increment for compatible new
features, and reserve `1.0.0` for the approved stable workflow.

## Output compatibility

New result folders use `results/*.mat` and do not emit the historical
`voltage_results.mat`, `trace_results.mat`, or `dual_results.mat` containers.
AAA can read old result folders through the single compatibility loader;
historical scripts that only understand the old filenames cannot read new
AAA output without that loader. Numerical algorithms and trace-stage
provenance remain unchanged.

Before approving the migration boundary, read `ARCHITECTURE.md`,
`IMPLEMENTATION_DECISIONS.md`, and `MIGRATION_STATUS.md`. The last file
distinguishes the completed App/control foundation from scientific sections
that remain explicitly excluded or require real-data parity fixtures.

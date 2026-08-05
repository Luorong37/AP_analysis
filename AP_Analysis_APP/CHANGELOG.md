# AAA change log

This file records user-visible releases of AP Analysis App (AAA). Version
numbers follow `major.minor.patch`:

- `major`: incompatible request/result contract changes;
- `minor`: backward-compatible workflows, sections, or UI capabilities;
- `patch`: fixes that do not intentionally change numerical behavior.

## 0.1.0 - 2026-08-05

Initial versioned development release.

- Added the unified AP/Dual request, validation, dispatch, and manifest layer.
- Added shared Single/Record section execution and Record-average control.
- Added unified `AAA.result.v1` writing plus centralized legacy-result reading.
- Added bounded reuse discovery, preflight receipts, analysis logs, and App status.
- Added source provenance to every new manifest: App version, repository commit
  when readable, MATLAB release, dependency version, and a SHA-256 digest of
  the production source bundle.

Known limitation: AP event/FWHM/statistics parity is not yet released. This is
therefore a `0.x` development release rather than a stable `1.0.0` release.


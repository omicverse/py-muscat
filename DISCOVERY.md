# Discovery: muscat

Checked omicverse organization repository inventory on 2026-11-11: no standalone target port found.
Upstream commit: `861683175a1804430ced420397e72fff16fe92ed`.
Scope: aggregate_data, pb_ds; other package functions are excluded.

Reuse: NumPy/SciPy/pandas for matrix and SVD/NNLS operations. muscat uses existing pyedger and python-limma rather than rewriting either statistical engine.
RUVg already exists in omicverse PR #1035 and will be extracted with its reference tests rather than independently duplicated.
Nebulosa and scran pooling have external candidates and are separate reuse audits, not new ports.

Fixed numerical gates are in manifest.yaml, committed before implementation. Synthetic and real public fixtures will be distinguished. No biological accuracy claim follows from fixture-level equivalence.

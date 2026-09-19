# Minimal publication changes

- Only canonical workflow code, necessary configuration, concise documentation
  and available analytical inputs are tracked. Historical result archives,
  rendered figures, draws and fitted-network objects are excluded.
- Identical functional matrices and current metadata are stored once. An explicit
  staging utility supplies exact copies at runtime where legacy scripts expect them.
- The 809 runner reads sample/node order from its fixed CLR matrix. Its edge order
  is `b, a = numpy.tril_indices(n_nodes, -1)`, matching the archived column-wise
  upper-triangle traversal; ordinary row-wise `triu_indices` would be incorrect.
- The 99,999 profile-restricted permutations are regenerated with the unchanged
  schedule function and original NumPy seed 1032. The ridge, LIONESS, metric,
  resampling and inference functions are unchanged. Equality against the archived
  schedule and network weights is checked before publication.
- The archived-network equality assertion inside 809 is removed because the
  archived fitted weights are no longer required inputs. Independent numerical
  acceptance checks remain. `network_statistics.py` retains the imported
  inference functions from legacy803 and omits its unused historical CLI.
- Figure-package verification checks locally generated figure-input hashes;
  the old whole-output-package manifest is not part of this minimal repository.
- The input collector requests only the missing phyloseq/protein inputs and
  required metadata/QC controls. Unknown inventory sizes are allowed, but actual
  copied bytes and SHA-256 are always recorded. Original files are not modified.

R execution and a complete fresh gene-to-figure workflow were not repeated here.

## Numerical equivalence check

The three constructions reproduced all 1,014,900 sample-specific weights each
and both sample metrics with maximum absolute differences of 0.0 against the
recovered reference objects. All 99,999 permutation rows were identical. These
checks cover the full-data/leave-one-out calculations and permutation schedule;
the 499-draw full reconstruction inference was not repeated.

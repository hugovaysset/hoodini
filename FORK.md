# This fork

Forked from [pentamorfico/hoodini](https://github.com/pentamorfico/hoodini) at
`c58abd7`, to make the pipeline survive a large input. Upstream's defaults are
tuned for the tens-to-hundreds of targets the tool is normally pointed at; at
43,739 targets / 614,037 neighbour proteins several of them stop being slow and
start being impossible.

Every change is size-conditional. **A small run behaves exactly as it did
upstream** — same clustering, same HTML, same everything — which is what makes
these upstreamable rather than a private variant.

## What changed, and the measurement behind it

### Clustering parameters scale with the input

`pipeline/helpers/mmseqs_clustering.py` hardcoded `sensitivity=15`,
`cluster_steps=9` and five rounds of iterative profile search. MMseqs's own
help documents `-s` as "1.0 faster; 4.0 fast; 7.5 sensitive [4.000]" and
`--cluster-steps` as "from 1 to -s", so 15 and 9 are both off the end of the
scale. Measured on 614,037 neighbour proteins, 32 threads:

| clustering | wall | clusters |
|---|---|---|
| `-s 15 --cluster-steps 9` + 5 profile rounds | **17 h 33 m** | 128,931 |
| `-s 7.5 --cluster-steps 3` | 4 m 43 s | 189,994 |
| `-s 4 --cluster-steps 3` | 1 m 00 s | 212,886 |
| `linclust` | 7 s | 557,644 |

The profile rounds earn their keep — they merge 190k families down to 129k, so
a third of the remote homology only shows up that way — but at 223x the cost.
So `plan_for(n_seqs)` picks: the original deep settings below 50,000 sequences,
`-s 7.5 --cluster-steps 3` above that, `linclust` above 2,000,000. Every
parameter is now an argument, and passing one explicitly overrides the plan.

### Clustering can be resumed

`cluster_with_mmseqs(resume=True)` returns immediately if its output is already
there, and `--resume` stops `initialize` from deleting the output folder.
Upstream's `--force` does `shutil.rmtree(output_folder)` and there is no
checkpoint anywhere, so a crash in any later stage costs the whole clustering —
which is the difference between a large run being retryable and not.

### The standalone HTML is skippable

`write_data` base64s every parquet table into one HTML document. That is what
makes it portable and what stops it existing above a certain size: a
44-neighbourhood run gives 3.5 MB, and a 43,739-neighbourhood one projects to
roughly a quarter of a gigabyte, which no browser opens. `--no-html` writes
only the tables; `--html-max-mb` (default 64) skips the document rather than
spending ten minutes producing one that cannot be used. The size is checked
*before* anything is encoded.

### DefenseFinder is threaded, and bounded

`extra_tools/defensefinder.py` ran one single-process `defense-finder` over
every neighbour protein at once. It now passes `-w <threads>`, and refuses
above 200,000 proteins with a message rather than never returning.

### Two smaller fixes

- `_run_command` took a shell string and called `sys.exit(1)` on failure. It now
  takes a list — no quoting bug on a path with a space — and raises
  `MmseqsError`, because exiting the interpreter from inside a library takes the
  caller's error handling away from it. MMseqs also finally receives
  `--threads`; `--num-threads` never reached it before.
- A TOML config of top-level scalars crashed with `'bool' object has no
  attribute 'items'`, naming neither the file nor the line: `cli.run` assumed
  every top-level key was a table, while `load_default_config` already accepted
  both shapes. They now agree.

## Running it without installing

The Aleph conda environment is root-owned, so `pip install -e` is not
available. `bin/hoodini` puts this checkout on `PYTHONPATH` ahead of
site-packages and runs the env's own interpreter, so the fork shadows the
installed package while polars, ete3, jinja2 and the rest still come from the
environment that has them:

```bash
export HOODINI_BIN=/home/ubuntu/hoodini/bin/hoodini
aleph_hoodini --inputsheet sheet.tsv --output out --config big.toml
```

`aleph_hoodini` forwards `--config`, which is how `--no-html` and `--resume`
reach Hoodini through the wrapper.

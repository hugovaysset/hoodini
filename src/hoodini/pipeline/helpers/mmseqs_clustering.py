"""Cluster neighbour proteins with MMseqs2.

Scaling
-------
The original defaults here — ``sensitivity=15``, ``cluster_steps=9`` and five
rounds of iterative profile search — are tuned for the few hundred proteins a
sixty-target neighbourhood set produces, where the whole thing costs seconds
and the extra remote homology is free. They do not survive a large input.

Measured on 614,037 neighbour proteins (43,739 SIR2 targets), 32 threads:

===========================================  =========  ============
clustering                                   wall       clusters
===========================================  =========  ============
``-s 15 --cluster-steps 9`` + 5 profile       17h 33m    128,931
``-s 7.5 --cluster-steps 3``                  4m 43s     189,994
``-s 4 --cluster-steps 3``                    1m 00s     212,886
``linclust``                                  7s         557,644
===========================================  =========  ============

The profile rounds are doing real work — they merge 190k families down to 129k,
so a third of the remote homology only shows up that way — but at 223 times the
cost. So the parameters become arguments with **size-aware defaults** rather
than constants: small inputs keep exactly the behaviour they had, and large
ones get a setting that finishes.

Note also that MMseqs's own help caps the scale it documents at 7.5
("Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [4.000]") and describes
``--cluster-steps`` as "from 1 to -s". ``-s 15`` is off the end of both.
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import polars as pl

from hoodini.utils.logging_utils import error, info

#: Above this many sequences the iterative profile search is skipped. It is
#: quadratic-ish in practice and it is what turns minutes into a day.
PROFILE_MAX_SEQS = 50_000

#: Above this, even a cascaded all-vs-all is too slow and clustering falls back
#: to linclust, which is linear in the number of sequences.
LINCLUST_MIN_SEQS = 2_000_000


class MmseqsError(RuntimeError):
    """An MMseqs2 command failed."""


def _run_command(command, log=None):
    """Run one MMseqs2 command.

    A list, not a string: ``shell=True`` on paths that can contain a space is a
    quoting bug waiting to happen. And an exception rather than ``sys.exit(1)``
    — this is a library, and exiting the interpreter from inside one takes the
    caller's own error handling away from it.
    """
    if isinstance(command, str):
        command = command.split()
    command = [str(c) for c in command]
    if log is not None:
        with open(log, "a") as fh:
            fh.write("$ " + " ".join(command) + "\n")
            result = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT, text=True)
        if result.returncode != 0:
            raise MmseqsError(
                f"mmseqs failed (exit {result.returncode}): {' '.join(command)}\n"
                f"  see {log}"
            )
        return ""
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        tail = "\n  ".join((result.stderr or "").strip().splitlines()[-6:])
        raise MmseqsError(
            f"mmseqs failed (exit {result.returncode}): {' '.join(command)}\n  {tail}"
        )
    return result.stdout


def read_output_to_df(output_file):
    df = pl.read_csv(output_file, separator="\t", header=None)
    return df


def count_sequences(fasta: Path) -> int:
    """How many records are in a FASTA, without holding it in memory."""
    n = 0
    with open(fasta, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                n += 1
    return n


def plan_for(n_seqs: int) -> dict:
    """Clustering parameters appropriate to an input of this size.

    Explicit rather than a smooth curve, because the three regimes really are
    different algorithms and it should be obvious from the log which one ran.
    """
    if n_seqs >= LINCLUST_MIN_SEQS:
        return {"mode": "linclust", "sensitivity": None, "cluster_steps": None, "max_steps": 0}
    if n_seqs > PROFILE_MAX_SEQS:
        return {"mode": "cluster", "sensitivity": 7.5, "cluster_steps": 3, "max_steps": 0}
    return {"mode": "cluster", "sensitivity": 15, "cluster_steps": 9, "max_steps": 5}


def cluster_with_mmseqs(
    fasta,
    temp_folder,
    max_steps=None,
    sensitivity=None,
    cluster_mode=1,
    cluster_steps=None,
    cov_mode=0,
    coverage=0.7,
    output=None,
    threads=None,
    log=None,
    resume=True,
):
    """Cluster ``fasta`` into families, writing a two-column TSV to ``output``.

    ``max_steps``, ``sensitivity`` and ``cluster_steps`` default to `plan_for`
    on the input's size. Passing any of them explicitly overrides that, so an
    existing caller that set them keeps its behaviour exactly.

    ``resume`` returns immediately if ``output`` already exists and is
    non-empty. Clustering is by far the longest stage, it is deterministic, and
    losing hours of it to a crash in a later stage is the difference between a
    pipeline that can be run on a large input and one that cannot.
    """
    fasta = Path(fasta)
    temp_folder = Path(temp_folder)
    output = Path(output) if output is not None else None

    if resume and output is not None and output.exists() and output.stat().st_size > 0:
        info(f"↩️\tReusing existing clustering at {output}")
        return

    n_seqs = count_sequences(fasta)
    plan = plan_for(n_seqs)
    # An explicit sensitivity is a request for a real cascaded cluster run, so
    # it overrides the linclust regime as well as the number it names.
    asked_for_cluster = sensitivity is not None
    if sensitivity is None:
        sensitivity = plan["sensitivity"]
    if cluster_steps is None:
        cluster_steps = plan["cluster_steps"] or 3
    if max_steps is None:
        max_steps = plan["max_steps"]
    mode = "cluster" if (asked_for_cluster or plan["mode"] == "cluster") else "linclust"

    if threads is None:
        threads = os.cpu_count() or 1
    thread_args = ["--threads", str(int(threads))]

    info(
        f"🧬\tClustering {n_seqs:,} proteins with mmseqs {mode}"
        + (f" -s {sensitivity} --cluster-steps {cluster_steps}" if mode == "cluster" else "")
        + (f", {max_steps} profile round(s)" if max_steps else "")
        + f", {threads} threads"
    )

    temp_folder.mkdir(parents=True, exist_ok=True)
    dbname = "temp_db"
    _run_command(["mmseqs", "createdb", fasta, temp_folder / dbname], log)

    if mode == "linclust":
        _run_command(
            ["mmseqs", "linclust", temp_folder / dbname, temp_folder / (dbname + "_clu"),
             temp_folder, "--cluster-mode", cluster_mode, "--cov-mode", cov_mode,
             "-c", coverage, *thread_args],
            log,
        )
    else:
        _run_command(
            ["mmseqs", "cluster", temp_folder / dbname, temp_folder / (dbname + "_clu"),
             temp_folder, "-s", sensitivity, "--cluster-mode", cluster_mode,
             "--cluster-steps", cluster_steps, "--cov-mode", cov_mode, "-c", coverage,
             *thread_args],
            log,
        )

    if max_steps:
        _run_command(
            ["mmseqs", "createsubdb", temp_folder / (dbname + "_clu"), temp_folder / dbname,
             temp_folder / (dbname + "_clu_repseq")],
            log,
        )
        input_db = f"{dbname}_clu_repseq"
        for step in range(1, int(max_steps) + 1):
            _run_command(
                ["mmseqs", "search", temp_folder / input_db, temp_folder / input_db,
                 temp_folder / f"search_step_{step}", temp_folder, "--add-self-matches",
                 *thread_args],
                log,
            )
            _run_command(
                ["mmseqs", "result2profile", temp_folder / (dbname + "_clu_repseq"),
                 temp_folder / (dbname + "_clu_repseq"), temp_folder / f"search_step_{step}",
                 temp_folder / f"search_step_{step}_profile", *thread_args],
                log,
            )
            _run_command(
                ["mmseqs", "search", temp_folder / f"search_step_{step}_profile",
                 temp_folder / input_db, temp_folder / f"search_step_{step}_pp", temp_folder,
                 "--add-self-matches", *thread_args],
                log,
            )
            _run_command(
                ["mmseqs", "clust", temp_folder / f"search_step_{step}_profile",
                 temp_folder / f"search_step_{step}_pp",
                 temp_folder / f"search_step_{step}_pp_clu", *thread_args],
                log,
            )
            _run_command(
                ["mmseqs", "createsubdb", temp_folder / f"search_step_{step}_pp_clu",
                 temp_folder / dbname, temp_folder / f"search_step_{step}_pp_clu_repseq"],
                log,
            )
            input_db = f"search_step_{step}_pp_clu_repseq"

        cluster_files = [
            temp_folder / f"search_step_{step}_pp_clu" for step in range(1, int(max_steps) + 1)
        ]
        _run_command(
            ["mmseqs", "mergeclusters", temp_folder / dbname, temp_folder / "deep_cluster_db",
             temp_folder / (dbname + "_clu"), *cluster_files],
            log,
        )
        final_db = temp_folder / "deep_cluster_db"
    else:
        final_db = temp_folder / (dbname + "_clu")

    _run_command(
        ["mmseqs", "createtsv", temp_folder / dbname, temp_folder / dbname, final_db, output,
         *thread_args],
        log,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Cluster protein sequences using MMseqs2 with optional profile rounds."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input FASTA file containing protein sequences"
    )
    parser.add_argument(
        "-t", "--temp-folder", required=True, help="Folder to save intermediary files"
    )
    parser.add_argument(
        "--max-steps",
        type=int,
        default=None,
        help="Iterative profile-search rounds (default: by input size)",
    )
    parser.add_argument(
        "-s", "--sensitivity", type=float, default=None,
        help="Sensitivity for clustering (default: by input size)"
    )
    parser.add_argument("--cluster-mode", type=int, default=1, help="Clustering mode (default: 1)")
    parser.add_argument(
        "--cluster-steps", type=int, default=None,
        help="Cascaded clustering steps (default: by input size)"
    )
    parser.add_argument(
        "--cov-mode", type=int, default=0, help="Coverage mode for clustering (default: 0)"
    )
    parser.add_argument(
        "-c", "--coverage", type=float, default=0.7, help="Coverage for clustering (default: 0.7)"
    )
    parser.add_argument("--threads", type=int, default=None, help="Threads (default: all cores)")
    parser.add_argument(
        "--no-resume", action="store_true", help="Recluster even if the output already exists"
    )
    parser.add_argument("-o", "--output", required=True, help="Output file name")
    parser.add_argument(
        "--keep-temp", action="store_true", help="Keep the MMseqs2 temp folder"
    )
    args = parser.parse_args()

    try:
        cluster_with_mmseqs(
            args.input,
            args.temp_folder,
            args.max_steps,
            args.sensitivity,
            args.cluster_mode,
            args.cluster_steps,
            args.cov_mode,
            args.coverage,
            args.output,
            threads=args.threads,
            resume=not args.no_resume,
        )
    except MmseqsError as exc:
        error(str(exc))
        sys.exit(1)
    if not args.keep_temp:
        shutil.rmtree(args.temp_folder, ignore_errors=True)


if __name__ == "__main__":
    main()

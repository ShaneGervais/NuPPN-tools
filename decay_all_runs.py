#!/usr/bin/env python3
"""Decay every PPN run in a nova case and mirror runs/ into decay_runs/."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


DEFAULT_DECAY_TIME_SECONDS = 2.0 * 60.0 * 60.0
TEXT_INPUTS = {
    "Makefile",
    "isotopedatabase.txt",
    "isotopedatabase_all.txt",
    "isotopedatabase_cf.txt",
    "networksetup.txt",
    "ppn_frame.input",
    "ppn_physics.input",
    "ppn_solver.input",
}


@dataclass(frozen=True)
class DecayJob:
    source_run: Path
    decay_run: Path
    source_iso_massf: Path
    final_cycle: int


@dataclass(frozen=True)
class DecayResult:
    job: DecayJob
    returncode: int
    output: Path
    log: Path
    status: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Mirror nova_case/runs into nova_case/decay_runs and run PPN decay "
            "for each run's final iso_massf file."
        )
    )
    parser.add_argument("nova_case", type=Path, help="Nova case directory containing runs/")
    parser.add_argument(
        "--runs-dir",
        type=Path,
        default=None,
        help="Source runs directory. Default: nova_case/runs",
    )
    parser.add_argument(
        "--decay-runs-dir",
        type=Path,
        default=None,
        help="Output decay runs directory. Default: nova_case/decay_runs",
    )
    parser.add_argument(
        "--decay-time",
        type=float,
        default=DEFAULT_DECAY_TIME_SECONDS,
        help="Decay time in seconds. Default: 7200 seconds, or 2 hours.",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=4,
        help="Number of ppn.exe decay runs to execute in parallel. Default: 4",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the jobs that would be built and executed without writing files.",
    )
    parser.add_argument(
        "--no-run",
        action="store_true",
        help="Build decay run directories and manifest without executing ppn.exe.",
    )
    return parser.parse_args()


def cycle_name(cycle: int) -> str:
    return f"{cycle:05d}"


def final_cycle_from_xtime(run_dir: Path) -> int:
    xtime = run_dir / "x-time.dat"
    if not xtime.is_file():
        raise FileNotFoundError(f"missing x-time.dat: {xtime}")

    final_cycle: int | None = None
    with xtime.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            final_cycle = int(stripped.split()[0])

    if final_cycle is None:
        raise ValueError(f"no cycle rows found in {xtime}")
    return final_cycle


def discover_source_runs(runs_dir: Path) -> list[Path]:
    return sorted(path.parent for path in runs_dir.rglob("ppn.exe") if path.is_file())


def resolve_relative_to(base: Path, path: Path) -> Path:
    path = path.expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def parse_iso_massf_rows(path: Path) -> list[dict[str, object]]:
    if not path.is_file():
        raise FileNotFoundError(f"missing iso_massf input: {path}")

    rows: list[dict[str, object]] = []
    with path.open() as handle:
        for line in handle:
            parts = line.split()
            if not parts or not parts[0].isdigit():
                continue
            rows.append(
                {
                    "z": int(float(parts[1])),
                    "a": int(float(parts[2])),
                    "abundance": float(parts[4]),
                    "label_parts": parts[5:],
                }
            )

    if not rows:
        raise ValueError(f"no abundance rows found in {path}")
    return rows


def write_post_abundance(source_iso_massf: Path, out_path: Path) -> None:
    rows = parse_iso_massf_rows(source_iso_massf)
    with out_path.open("w") as handle:
        for row in rows:
            z = int(row["z"])
            a = int(row["a"])
            abundance = float(row["abundance"])
            label_parts = list(row["label_parts"])

            if label_parts == ["NEUT"]:
                iso = "NEUT"
            elif label_parts == ["PROT"]:
                iso = "PROT"
            else:
                if len(label_parts) >= 2:
                    element = label_parts[0].lower()
                else:
                    match = re.match(r"^([A-Za-z]+)(\d+)$", label_parts[0])
                    if not match:
                        raise ValueError(
                            f"cannot convert isotope label from {source_iso_massf}: {label_parts}"
                        )
                    element = match.group(1).lower()
                iso = f"{element:<2s}{a:3d}"

            handle.write(f"{z:3d} {iso:<5s}         {abundance:16.10E}\n")


def copy_or_link(src: Path, dst: Path) -> None:
    if src.is_symlink():
        target = src.resolve()
        dst.symlink_to(Path(os.path.relpath(target, start=dst.parent)))
    else:
        shutil.copy2(src, dst)


def ensure_symlink(link_path: Path, target: Path) -> None:
    if link_path.exists() or link_path.is_symlink():
        return
    link_path.symlink_to(Path(os.path.relpath(target.resolve(), start=link_path.parent)))


def copy_decay_inputs(source_run: Path, decay_run: Path) -> None:
    for name in sorted(TEXT_INPUTS):
        src = source_run / name
        if src.exists():
            copy_or_link(src, decay_run / name)

    for src in sorted(source_run.glob("*.input")):
        dst = decay_run / src.name
        if not dst.exists():
            copy_or_link(src, dst)

    ppn_exe = source_run / "ppn.exe"
    if not ppn_exe.exists():
        raise FileNotFoundError(f"missing ppn.exe in source run: {ppn_exe}")
    (decay_run / "ppn.exe").symlink_to(Path(os.path.relpath(ppn_exe.resolve(), start=decay_run)))

    npdata = source_run.parent / "NPDATA"
    if npdata.exists():
        ensure_symlink(decay_run.parent / "NPDATA", npdata)


def update_namelist(text: str, replacements: dict[str, str]) -> str:
    lines = text.splitlines()
    found: set[str] = set()
    output: list[str] = []
    replacement_keys = {key.lower(): value for key, value in replacements.items()}

    for line in lines:
        stripped = line.strip()
        if stripped == "/":
            for key, value in replacements.items():
                if key.lower() not in found:
                    output.append(f"        {key} = {value}")
            output.append(line)
            continue

        match = re.match(r"^(\s*)([A-Za-z_][A-Za-z0-9_]*)\s*=", line)
        if match and match.group(2).lower() in replacement_keys:
            key = match.group(2)
            found.add(key.lower())
            output.append(f"{match.group(1)}{key} = {replacement_keys[key.lower()]}")
        else:
            output.append(line)

    return "\n".join(output) + "\n"


def fortran_float(value: float) -> str:
    return f"{value:.10E}".replace("E", "d")


def patch_decay_inputs(decay_run: Path, decay_time: float) -> None:
    frame = decay_run / "ppn_frame.input"
    physics = decay_run / "ppn_physics.input"
    if not frame.is_file():
        raise FileNotFoundError(f"missing ppn_frame.input in decay run: {frame}")
    if not physics.is_file():
        raise FileNotFoundError(f"missing ppn_physics.input in decay run: {physics}")

    frame.write_text(
        update_namelist(
            frame.read_text(),
            {
                "nsource": "0",
                "iabuini": "11",
                "ini_filename": "'post_abundance.DAT'",
                "iplot_flux_option": "0",
                "i_flux_integrated": "0",
            },
        )
    )
    physics.write_text(
        update_namelist(
            physics.read_text(),
            {
                "decay": ".true.",
                "decay_time": fortran_float(decay_time),
                "detailed_balance": ".false.",
            },
        )
    )


def build_jobs(runs_dir: Path, decay_runs_dir: Path, source_runs: list[Path]) -> list[DecayJob]:
    jobs: list[DecayJob] = []
    for source_run in source_runs:
        final_cycle = final_cycle_from_xtime(source_run)
        source_iso_massf = source_run / f"iso_massf{cycle_name(final_cycle)}.DAT"
        if not source_iso_massf.is_file():
            raise FileNotFoundError(f"missing final iso_massf file: {source_iso_massf}")
        relative = source_run.relative_to(runs_dir)
        jobs.append(
            DecayJob(
                source_run=source_run,
                decay_run=decay_runs_dir / relative,
                source_iso_massf=source_iso_massf,
                final_cycle=final_cycle,
            )
        )
    return jobs


def build_decay_run(nova_case: Path, job: DecayJob, decay_time: float) -> None:
    job.decay_run.mkdir(parents=True)
    copy_decay_inputs(job.source_run, job.decay_run)
    write_post_abundance(job.source_iso_massf, job.decay_run / "post_abundance.DAT")
    patch_decay_inputs(job.decay_run, decay_time)
    (job.decay_run / "decay_source.txt").write_text(
        "\n".join(
            [
                f"nova_case = {nova_case}",
                f"source_run = {job.source_run}",
                f"source_iso_massf = {job.source_iso_massf}",
                f"source_final_cycle = {job.final_cycle}",
                f"decay_time_seconds = {decay_time:.10g}",
                f"generated = {datetime.now().isoformat(timespec='seconds')}",
                "",
            ]
        )
    )


def ppn_env() -> dict[str, str]:
    env = os.environ.copy()
    for key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        env[key] = "1"
    return env


def run_ppn_decay(job: DecayJob, env: dict[str, str]) -> DecayResult:
    log_path = job.decay_run / "decay.log"
    with log_path.open("w") as log:
        completed = subprocess.run(
            ["./ppn.exe"],
            cwd=job.decay_run,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
        )

    output = job.decay_run / "iso_massfdecay.DAT"
    if completed.returncode == 0 and output.is_file():
        status = "ok"
    elif completed.returncode == 0:
        status = "missing_output"
    else:
        status = "failed"
    return DecayResult(job, completed.returncode, output, log_path, status)


def write_manifest(path: Path, jobs: list[DecayJob], results: list[DecayResult] | None = None) -> None:
    result_by_run = {} if results is None else {result.job.decay_run: result for result in results}
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "source_run",
                "decay_run",
                "source_iso_massf",
                "final_cycle",
                "output",
                "log",
                "returncode",
                "status",
            ],
        )
        writer.writeheader()
        for job in jobs:
            result = result_by_run.get(job.decay_run)
            writer.writerow(
                {
                    "source_run": job.source_run,
                    "decay_run": job.decay_run,
                    "source_iso_massf": job.source_iso_massf,
                    "final_cycle": job.final_cycle,
                    "output": "" if result is None else result.output,
                    "log": "" if result is None else result.log,
                    "returncode": "" if result is None else result.returncode,
                    "status": "built" if result is None else result.status,
                }
            )


def main() -> None:
    args = parse_args()
    nova_case = args.nova_case.resolve()
    if not nova_case.is_dir():
        raise SystemExit(f"nova_case does not exist: {nova_case}")

    runs_dir = (
        resolve_relative_to(nova_case, args.runs_dir)
        if args.runs_dir is not None
        else (nova_case / "runs").resolve()
    )
    decay_runs_dir = (
        resolve_relative_to(nova_case, args.decay_runs_dir)
        if args.decay_runs_dir is not None
        else (nova_case / "decay_runs").resolve()
    )
    if not runs_dir.is_dir():
        raise SystemExit(f"runs directory does not exist: {runs_dir}")
    if args.jobs < 1:
        raise SystemExit("--jobs must be at least 1")

    source_runs = discover_source_runs(runs_dir)
    if not source_runs:
        raise SystemExit(f"no ppn.exe files found below {runs_dir}")
    jobs = build_jobs(runs_dir, decay_runs_dir, source_runs)

    print(f"found {len(jobs)} source runs")
    print(f"decay time: {args.decay_time:.10g} seconds")
    print(f"runs dir: {runs_dir}")
    print(f"decay runs dir: {decay_runs_dir}")

    if args.dry_run:
        for job in jobs:
            print(f"{job.source_run} -> {job.decay_run}")
        return

    if decay_runs_dir.exists():
        print(f"WARNING: removing existing decay_runs directory: {decay_runs_dir}")
        shutil.rmtree(decay_runs_dir)

    decay_runs_dir.mkdir(parents=True)
    root_npdata = runs_dir / "NPDATA"
    if root_npdata.exists():
        ensure_symlink(decay_runs_dir / "NPDATA", root_npdata)

    for job in jobs:
        build_decay_run(nova_case, job, args.decay_time)

    manifest = decay_runs_dir / "decay_manifest.csv"
    write_manifest(manifest, jobs)
    print(f"built {len(jobs)} decay run directories")

    if args.no_run:
        print(f"wrote manifest: {manifest}")
        print("--no-run requested; not executing ppn.exe")
        return

    env = ppn_env()
    results: list[DecayResult] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        future_to_job = {executor.submit(run_ppn_decay, job, env): job for job in jobs}
        for idx, future in enumerate(as_completed(future_to_job), start=1):
            result = future.result()
            results.append(result)
            rel = result.job.decay_run.relative_to(decay_runs_dir)
            print(f"[{idx}/{len(jobs)}] {result.status}: {rel}")

    results.sort(key=lambda result: str(result.job.decay_run))
    write_manifest(manifest, jobs, results)
    failures = [result for result in results if result.status != "ok"]
    print(f"wrote manifest: {manifest}")
    print(f"successful decays: {len(results) - len(failures)}")
    print(f"failed decays: {len(failures)}")
    if failures:
        for result in failures[:20]:
            print(f"  {result.status}: {result.job.decay_run} log={result.log}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()

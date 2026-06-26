#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]

VALID_SOURCES = {
    "BASEL", "JINAR", "JINAC", "JINAV", "ILI01", "NACRL", "NACRR", "NACRU",
    "VITAL", "RVRSE", "KADON", "NETB1", "ODA94", "LMP00", "FFW85", "JBJ16",
    "NKK04", "MPG06", "NTRNO", "OOOOO",
}

ELEMENT_SYMBOLS = {
    "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE",
    "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", "K", "CA",
    "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
    "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR",
    "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN",
    "SB", "TE", "I", "XE", "CS", "BA", "LA", "CE", "PR", "ND",
    "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB",
    "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG",
    "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH",
    "PA", "U",
}

NETWORK_RE = re.compile(
    r"^\s*(\d+)\s+([TF])\s+(\d+)\s+(.{5})\s+\+\s+(\d+)\s+(.{5})"
    r"\s+->\s+(\d+)\s+(.{5})\s+\+\s+(\d+)\s+(.{5})\s+\S+\s+(\S+)\s+(\S+)\s+(\d+)"
)


@dataclass(frozen=True)
class Row:
    index: int
    active: bool
    reactant: tuple[int, str] | None
    projectile: tuple[int, str] | None
    product_1: tuple[int, str] | None
    product_2: tuple[int, str] | None
    source: str
    rtype: str
    line_no: int
    line: str


def nova_dir(name: str) -> Path:
    return PROJECT_ROOT / "nova_cases" / name


def parse_json_file(path: str | Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def species_id(text: str) -> tuple[int, str] | None:
    text = text.strip()
    if text == "PROT":
        return (1, "H")
    if text == "NEUT":
        return (1, "N")
    if text == "OOOOO":
        return (0, "G")
    match = re.match(r"([A-Z]+)\s*(\d+)$", text)
    if not match:
        return None
    return (int(match.group(2)), match.group(1))


def species_from_name(text: str) -> tuple[int, str]:
    match = re.match(r"(\d+)([A-Za-z]+)$", text)
    if not match:
        raise ValueError(f"Cannot parse isotope name {text!r}")
    symbol = match.group(2).upper()
    if symbol not in ELEMENT_SYMBOLS and len(symbol) > 1 and symbol[-1] in {"G", "M"}:
        base = symbol[:-1]
        if base in ELEMENT_SYMBOLS:
            symbol = base
    return (int(match.group(1)), symbol)


def reaction_from_name(name: str):
    parts = name.split("_")
    if len(parts) != 3:
        raise ValueError(f"Reaction name must look like '20Ne_pg_21Na': {name}")
    target, channel, product = parts
    channels = {
        "pg": ("(p,g)", (1, "H")),
        "pa": ("(p,a)", (1, "H")),
        "ag": ("(a,g)", (4, "HE")),
        "an": ("(a,n)", (4, "HE")),
        "ap": ("(a,p)", (4, "HE")),
    }
    if channel not in channels:
        raise ValueError(f"Unsupported reaction channel {channel!r} in {name}")
    rtype, projectile = channels[channel]
    return species_from_name(target), projectile, rtype, species_from_name(product)


def parse_networksetup(path: str | Path) -> list[Row]:
    rows: list[Row] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle):
            match = NETWORK_RE.match(line)
            if not match:
                continue
            c = match.groups()
            rows.append(Row(
                index=int(c[0]),
                active=c[1] == "T",
                reactant=species_id(c[3]),
                projectile=species_id(c[5]),
                product_1=species_id(c[7]),
                product_2=species_id(c[9]),
                source=c[10],
                rtype=c[11],
                line_no=line_no,
                line=line.strip(),
            ))
    return rows


def row_by_index(network: list[Row], index: int) -> Row | None:
    return next((row for row in network if row.index == index), None)


def matching_rows_for_reaction(reaction: dict[str, Any], network: list[Row]):
    name = reaction.get("network_name", reaction["name"])
    target, projectile, rtype, product = reaction_from_name(name)
    candidates = [
        row for row in network
        if row.reactant == target and row.projectile == projectile and row.rtype == rtype
        and (row.product_1 == product or row.product_2 == product)
    ]
    remapped = False
    if not candidates:
        candidates = [
            row for row in network
            if row.reactant == target and row.projectile == projectile and row.rtype == rtype
        ]
        remapped = bool(candidates)
    }
    remapped = remapped or bool(reaction.get("product_was_remapped", False))
    return candidates, remapped


def select_configured_row(reaction: dict[str, Any], candidates: list[Row]) -> Row | None:
    if "index" not in reaction:
        return None
    index = int(reaction["index"])
    return next((row for row in candidates if row.index == index), None)


def validate_reverse_index(reaction: dict[str, Any], network: list[Row], forward_index: int) -> int | None:
    if "reverse_index" not in reaction:
        return None
    reverse_index = int(reaction["reverse_index"])
    forward = row_by_index(network, forward_index)
    reverse = row_by_index(network, reverse_index)
    if forward is None:
        raise ValueError(f"{reaction['name']}: forward index {forward_index} was not found")
    if reverse is None:
        raise ValueError(f"{reaction['name']}: reverse_index {reverse_index} was not found")
    return reverse_index


def resolve_reaction_index(reaction: dict[str, Any], network: list[Row]) -> int:
    candidates, remapped = matching_rows_for_reaction(reaction, network)
    if not candidates:
        raise ValueError(f"{reaction['name']}: could not resolve reaction in networksetup.txt")
    selected = select_configured_row(reaction, candidates)
    if selected is None:
        raise ValueError(f"{reaction['name']}: configured index {reaction.get('index')} was not found")
    print(f"{reaction['name']}: using index {selected.index} ({selected.source})")
    if remapped:
        print(f"  note: product was remapped by network boundaries: {selected.line}")
    return selected.index


def copy_ppn(ppn_dir: str | Path, dest: str | Path) -> None:
    ppn_dir = Path(ppn_dir)
    dest = Path(dest)
    if dest.exists() or dest.is_symlink():
        if dest.is_dir() and not dest.is_symlink():
            shutil.rmtree(dest)
        else:
            dest.unlink()
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(ppn_dir, dest, symlinks=True)
    npdata_src = (ppn_dir / ".." / "NPDATA").resolve()
    for link in (dest.parent / "NPDATA", dest / "NPDATA"):
        if not link.exists() and not link.is_symlink():
            link.symlink_to(npdata_src, target_is_directory=True)


def set_starlib_option(path: str | Path, option: int) -> None:
    path = Path(path)
    lines = path.read_text().splitlines(keepends=True)
    for i, line in enumerate(lines):
        if re.match(r"(?i)^\s*starlib_option\s*=", line):
            lines[i] = re.sub(r"(?i)^(\s*starlib_option\s*=\s*)\d+", rf"\g<1>{option}", line, count=1)
            path.write_text("".join(lines))
            return
    for i, line in enumerate(lines):
        if line.strip() == "/":
            lines.insert(i, f"        starlib_option = {option}\n")
            path.write_text("".join(lines))
            return
    raise ValueError(f"Could not find ppn_physics namelist terminator in {path}")


def write_physics_input(ppn_dir: str | Path, run_dir: str | Path, indices: list[int], factors: list[float]) -> None:
    template = Path(ppn_dir) / "ppn_physics.input"
    output = Path(run_dir) / "ppn_physics.input"
    lines = template.read_text().splitlines(keepends=True)
    new_lines: list[str] = []
    for line in lines:
        if line.strip() == "/":
            for i, (idx, factor) in enumerate(zip(indices, factors), start=1):
                new_lines.append(f"        rate_index({i}) = {idx}\n")
                new_lines.append(f"        rate_factor({i}) = {factor:.16E}\n")
        new_lines.append(line)
    output.write_text("".join(new_lines))


def create_factored_runs(nova: str, baseline_only: bool = False, dry_run: bool = False, runs_name: str = "runs") -> Path:
    base = nova_dir(nova)
    ppn_dir = base / "ppn"
    runs_dir = base / runs_name
    config = parse_json_file(base / "config" / "reaction_plan.json")
    network = parse_networksetup(ppn_dir / "networksetup.txt")
    default_factors = config.get("default_factors", [])
    baseline_dir = runs_dir / "baseline"
    if dry_run:
        print(f"Would build baseline run in {baseline_dir}")
    else:
        runs_dir.mkdir(parents=True, exist_ok=True)
        copy_ppn(ppn_dir, baseline_dir)
        print(f"Built baseline run in {baseline_dir}")
    if baseline_only:
        return runs_dir

    created = 0
    for reaction in config["reactions"]:
        index = resolve_reaction_index(reaction, network)
        reverse_index = validate_reverse_index(reaction, network, index)
        for factor in reaction.get("factors", default_factors):
            run_dir = runs_dir / reaction["name"] / f"fact_{factor}"
            indices = [index]
            factors = [float(factor)]
            if reverse_index is not None:
                indices.append(reverse_index)
                factors.append(float(factor))
            if not dry_run:
                copy_ppn(ppn_dir, run_dir)
                write_physics_input(ppn_dir, run_dir, indices, factors)
            created += 1
    print(f"{'Would build' if dry_run else 'Built'} {created} factored runs in {runs_dir}")
    return runs_dir


def list_ppn_executables(runs_dir: str | Path, baseline_only: bool = False) -> list[Path]:
    runs_dir = Path(runs_dir)
    if baseline_only:
        exe = runs_dir / "baseline" / "ppn.exe"
        if not exe.is_file():
            raise FileNotFoundError(f"Missing baseline executable: {exe}")
        return [exe]
    exes = sorted(path for path in runs_dir.rglob("ppn.exe") if path.is_file())
    if not exes:
        raise FileNotFoundError(f"No ppn.exe files found under {runs_dir}")
    return exes


def run_one_ppn(exe: Path, runs_dir: Path, logs_dir: Path) -> None:
    run_dir = exe.parent
    name = str(run_dir.relative_to(runs_dir)).replace(os.sep, "_")
    logfile = logs_dir / f"{name}.log"
    logs_dir.mkdir(parents=True, exist_ok=True)
    start = time.time()
    with open(logfile, "w", encoding="utf-8") as log:
        log.write("===========================================\n")
        log.write(f"Running {run_dir}\n")
        log.write("===========================================\n")
        process = subprocess.run(["./ppn.exe"], cwd=run_dir, stdout=log, stderr=log)
    if process.returncode != 0:
        raise RuntimeError(f"Command failed in {run_dir}; see {logfile}")
    elapsed = round(time.time() - start)
    with open(logfile, "a", encoding="utf-8") as log:
        log.write(f"Finished {run_dir} in {elapsed} seconds\n")
    print(f"Finished {name} in {elapsed}s")


def run_parallel(exes: list[Path], runs_dir: Path, logs_dir: Path, jobs: int) -> None:
    from concurrent.futures import ThreadPoolExecutor, as_completed

    failures: list[tuple[Path, BaseException]] = []
    with ThreadPoolExecutor(max_workers=min(jobs, len(exes))) as pool:
        futures = {pool.submit(run_one_ppn, exe, runs_dir, logs_dir): exe for exe in exes}
        for future in as_completed(futures):
            exe = futures[future]
            try:
                future.result()
            except BaseException as exc:
                failures.append((exe, exc))
    if failures:
        print("Failures:")
        for exe, exc in failures:
            print(f"  {exe}\n    {exc}")
        raise RuntimeError(f"{len(failures)} ppn jobs failed")


def set_thread_env() -> None:
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ.setdefault(key, "1")


def final_iso_massf(run_dir: str | Path) -> Path | None:
    files = []
    for path in Path(run_dir).glob("iso_massf*.DAT"):
        match = re.match(r"iso_massf(\d+)\.DAT$", path.name)
        if match:
            files.append((int(match.group(1)), path))
    return max(files)[1] if files else None


def parse_iso_massf(path: str | Path) -> dict[str, float]:
    values: dict[str, float] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("H "):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                x = float(parts[4].replace("D", "E").replace("d", "e"))
            except ValueError:
                continue
            iso = " ".join(parts[5:])
            if iso == "PROT":
                iso = "H-1"
            elif iso != "NEUT":
                iso_parts = iso.split()
                if len(iso_parts) >= 2:
                    iso = f"{iso_parts[0].upper()}-{int(float(iso_parts[1]))}"
            values[iso] = x
    return values

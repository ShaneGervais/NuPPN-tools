#!/usr/bin/env python3
"""Extract Iliadis et al. 2001 Tables 3-9 and compare with NuPPN netgen.

The PDF text uses E[ for negative exponents and E] for positive exponents.
This script normalizes those markers and writes tab-separated .dat files.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import math
from pathlib import Path
import re
import subprocess
import tempfile


TABLE_COLUMNS = {
    3: [
        "20Ne(p,g)21Na",
        "21Ne(p,g)22Na",
        "22Ne(p,g)23Na",
        "20Na(p,g)21Mg",
        "21Na(p,g)22Mg",
        "22Na(p,g)23Mg",
        "23Na(p,g)24Mg",
        "23Na(p,a)20Ne",
    ],
    4: [
        "22Mg(p,g)23Al",
        "23Mg(p,g)24Al",
        "24Mg(p,g)25Al",
        "25Mg(p,g)26Al_total",
        "25Mg(p,g)26Alg",
        "25Mg(p,g)26Alm",
        "26Mg(p,g)27Al",
        "23Al(p,g)24Si",
    ],
    5: [
        "24Al(p,g)25Si",
        "25Al(p,g)26Si",
        "26Alg(p,g)27Si",
        "26Alm(p,g)27Si",
        "27Al(p,g)28Si",
        "27Al(p,a)24Mg",
        "26Si(p,g)27P",
        "27Si(p,g)28P",
    ],
    6: [
        "28Si(p,g)29P",
        "29Si(p,g)30P",
        "30Si(p,g)31P",
        "27P(p,g)28S",
        "28P(p,g)29S",
        "29P(p,g)30S",
        "30P(p,g)31S",
        "31P(p,g)32S",
    ],
    7: [
        "31P(p,a)28Si",
        "30S(p,g)31Cl",
        "31S(p,g)32Cl",
        "32S(p,g)33Cl",
        "33S(p,g)34Cl",
        "34S(p,g)35Cl",
        "31Cl(p,g)32Ar",
        "32Cl(p,g)33Ar",
    ],
    8: [
        "33Cl(p,g)34Ar",
        "34Cl(p,g)35Ar",
        "35Cl(p,g)36Ar",
        "35Cl(p,a)32S",
        "34Ar(p,g)35K",
        "35Ar(p,g)36K",
        "36Ar(p,g)37K",
        "35K(p,g)36Ca",
    ],
    9: [
        "36K(p,g)37Ca",
        "37K(p,g)38Ca",
        "38K(p,g)39Ca",
        "39K(p,g)40Ca",
        "39K(p,a)36Ar",
        "39Ca(p,g)40Sc",
        "40Ca(p,g)41Sc",
    ],
}

TABLE_CAPTION = "Recommended reaction rates N_A<sigma v>_gs in cm3 mol-1 s-1"
SOURCE_NOTE = "Iliadis et al. 2001, ApJS 134:151-171"
SYMBOLS = {
    "NE": "Ne",
    "NA": "Na",
    "MG": "Mg",
    "AL": "Al",
    "SI": "Si",
    "P": "P",
    "S": "S",
    "CL": "Cl",
    "AR": "Ar",
    "K": "K",
    "CA": "Ca",
    "SC": "Sc",
    "HE": "He",
}


@dataclass(frozen=True)
class ParsedTable:
    number: int
    columns: list[str]
    rows: list[tuple[float, list[str]]]


@dataclass(frozen=True)
class ArticlePoint:
    t9: float
    rate: float | None
    table: int
    column: str


@dataclass(frozen=True)
class NetgenBlock:
    reaction: str
    t9: list[float]
    rates: list[float]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_pdf_path() -> Path:
    return repo_root() / "ppn_nova" / "references" / "Iliadis-ApJS-2001.pdf"


def default_output_dir() -> Path:
    return repo_root() / "ppn_nova" / "references"


def default_netgen_path() -> Path:
    return repo_root().parent / "physics" / "NPDATA" / "netgen" / "netgen_iliadis2001_log_100.txt"


def extract_pdf_layout(pdf_path: Path) -> list[str]:
    with tempfile.NamedTemporaryFile(prefix="iliadis2001.", suffix=".layout.txt", delete=False) as handle:
        layout_path = Path(handle.name)

    subprocess.run(
        ["pdftotext", "-layout", str(pdf_path), str(layout_path)],
        check=True,
    )
    return layout_path.read_text(errors="replace").splitlines()


def normalize_rate_token(token: str) -> str:
    if token == "...":
        return token
    return token.replace("E[", "E-").replace("E]", "E+").replace("D", "E")


def parse_rate_token(token: str) -> float | None:
    normalized = normalize_rate_token(token)
    if normalized == "...":
        return None
    return float(normalized)


def table_bounds(lines: list[str], table_number: int) -> tuple[int, int]:
    table_re = re.compile(rf"\bTABLE\s+{table_number}\b")
    start = next(i for i, line in enumerate(lines) if table_re.search(line))
    end = next(
        i
        for i in range(start + 1, len(lines))
        if "NOTE" in lines[i] and f"Table {table_number}" in lines[i]
    )
    return start, end


def parse_tables(lines: list[str]) -> dict[int, ParsedTable]:
    parsed: dict[int, ParsedTable] = {}
    for table_number, columns in TABLE_COLUMNS.items():
        start, end = table_bounds(lines, table_number)
        rows: list[tuple[float, list[str]]] = []
        for line in lines[start:end]:
            tokens = [token for token in line.strip().split() if token != "."]
            if not tokens:
                continue
            try:
                t9 = float(tokens[0])
            except ValueError:
                continue

            values = tokens[1:]
            if len(values) != len(columns):
                continue
            rows.append((t9, [normalize_rate_token(value) for value in values]))

        if len(rows) != 31:
            raise ValueError(f"Table {table_number} parsed {len(rows)} rows, expected 31")
        parsed[table_number] = ParsedTable(table_number, columns, rows)
    return parsed


def write_article_tables(tables: dict[int, ParsedTable], output_dir: Path, pdf_path: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for table in tables.values():
        path = output_dir / f"iliadis2001_table{table.number:02d}.dat"
        with path.open("w") as handle:
            handle.write(f"# Source: {SOURCE_NOTE}; {pdf_path}\n")
            handle.write(f"# Caption: Table {table.number}. {TABLE_CAPTION}\n")
            handle.write("# Extracted from the local PDF with pdftotext -layout.\n")
            handle.write("# Reaction notation normalized: g=gamma, a=alpha.\n")
            handle.write("# Missing entries are kept as ... from the source table.\n")
            if table.number == 4:
                handle.write("# 25Mg(p,g)26Al_total/g/m follow the source t/g/m columns.\n")
            handle.write("\t".join(["T9_GK", *table.columns]) + "\n")
            for t9, values in table.rows:
                handle.write("\t".join([f"{t9:g}", *values]) + "\n")


def comparison_key(column: str) -> str:
    return column.replace("_total", "")


def article_series(tables: dict[int, ParsedTable]) -> dict[str, list[ArticlePoint]]:
    series: dict[str, list[ArticlePoint]] = {}
    for table in tables.values():
        for column_index, column in enumerate(table.columns):
            key = comparison_key(column)
            points = series.setdefault(key, [])
            for t9, values in table.rows:
                points.append(
                    ArticlePoint(
                        t9=t9,
                        rate=parse_rate_token(values[column_index]),
                        table=table.number,
                        column=column,
                    )
                )
    return series


def parse_species_line(line: str) -> tuple[bool, str | None]:
    match = re.match(r"#\s+(\d+)\s+([A-Z]+)\s*([0-9]+[gm]?)?\s*$", line)
    if match is None:
        return False, None

    count = int(match.group(1))
    symbol = match.group(2)
    mass = match.group(3)
    if count == 0 or symbol == "OOOOO":
        return True, None
    if symbol == "PROT":
        return True, "p"

    element = SYMBOLS.get(symbol, symbol.title())
    if mass is None:
        return True, element
    mass_match = re.match(r"(\d+)([gm]?)", mass)
    if mass_match is None:
        return True, element
    return True, f"{mass_match.group(1)}{element}{mass_match.group(2)}"


def parse_netgen(path: Path) -> dict[str, NetgenBlock]:
    blocks: dict[str, NetgenBlock] = {}
    current_species: list[str | None] = []
    current_reaction: str | None = None
    current_t9: list[float] = []
    current_rates: list[float] = []

    def finish_current() -> None:
        nonlocal current_reaction, current_t9, current_rates
        if current_reaction is not None:
            blocks[current_reaction] = NetgenBlock(current_reaction, current_t9, current_rates)
        current_reaction = None
        current_t9 = []
        current_rates = []

    for line in path.read_text().splitlines():
        matched_species, species = parse_species_line(line)
        if matched_species:
            current_species.append(species)
            continue

        if line.startswith("#Qrad"):
            finish_current()
            if len(current_species) < 4:
                raise ValueError(f"Could not parse reaction header before: {line}")
            target, projectile, product_1, product_2 = current_species[-4:]
            if target is None or projectile != "p" or product_2 is None:
                raise ValueError(f"Unexpected netgen reaction species: {current_species[-4:]}")

            if product_1 is None:
                reaction_type = "(p,g)"
            elif product_1 == "4He":
                reaction_type = "(p,a)"
            else:
                reaction_type = "(?)"
            current_reaction = f"{target}{reaction_type}{product_2}"
            continue

        data_match = re.match(r"\s*([0-9.]+)\s+([.0-9DEe+\-]+)", line)
        if data_match is not None and current_reaction is not None:
            current_t9.append(float(data_match.group(1)) / 10.0)
            current_rates.append(float(data_match.group(2).replace("D", "E")))

    finish_current()
    return blocks


def positive_article_points(points: list[ArticlePoint]) -> list[tuple[float, float]]:
    return sorted((point.t9, point.rate) for point in points if point.rate is not None and point.rate > 0.0)


def log_interpolate_pairs(points: list[tuple[float, float]], t9: float) -> float | None:
    if len(points) < 2 or t9 < points[0][0] or t9 > points[-1][0]:
        return None
    for (t1, r1), (t2, r2) in zip(points, points[1:]):
        if t1 <= t9 <= t2:
            if t9 == t1:
                return r1
            if t9 == t2:
                return r2
            weight = (math.log10(t9) - math.log10(t1)) / (math.log10(t2) - math.log10(t1))
            return 10.0 ** (math.log10(r1) + weight * (math.log10(r2) - math.log10(r1)))
    return None


def log_interpolate_netgen(block: NetgenBlock, t9: float) -> float | None:
    return log_interpolate_pairs(list(zip(block.t9, block.rates)), t9)


def fmt_float(value: float | None, digits: int = 6) -> str:
    if value is None:
        return "..."
    return f"{value:.{digits}E}"


def fmt_plain(value: float | None) -> str:
    if value is None:
        return "..."
    return f"{value:.8g}"


def write_comparisons(
    series: dict[str, list[ArticlePoint]],
    netgen: dict[str, NetgenBlock],
    output_dir: Path,
    factor_threshold: float,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    log_threshold = math.log10(factor_threshold)
    summary_rows: list[list[str]] = []
    runtime_mismatches: list[list[str]] = []
    stored_mismatches: list[list[str]] = []

    for reaction in sorted(series):
        points = series[reaction]
        finite_points = positive_article_points(points)
        block = netgen.get(reaction)
        first_point = points[0]
        n_article_finite = len(finite_points)

        if block is None:
            summary_rows.append(
                [
                    reaction,
                    str(first_point.table),
                    first_point.column,
                    "missing_in_netgen",
                    str(n_article_finite),
                    "0",
                    "0",
                    "...",
                    "...",
                    "...",
                    "...",
                    "...",
                    "0",
                    "0",
                    "...",
                    "...",
                    "...",
                    "...",
                    "...",
                ]
            )
            continue

        runtime_count = 0
        runtime_bad = 0
        runtime_worst: tuple[float, float, float, float, float] | None = None
        for point in points:
            if point.rate is None or point.rate <= 0.0:
                continue
            netgen_rate = log_interpolate_netgen(block, point.t9)
            if netgen_rate is None or netgen_rate <= 0.0:
                continue
            runtime_count += 1
            log_ratio = math.log10(netgen_rate / point.rate)
            ratio = netgen_rate / point.rate
            candidate = (abs(log_ratio), point.t9, point.rate, netgen_rate, log_ratio, ratio)
            if runtime_worst is None or candidate[0] > runtime_worst[0]:
                runtime_worst = candidate
            if abs(log_ratio) > log_threshold:
                runtime_bad += 1
                runtime_mismatches.append(
                    [
                        reaction,
                        str(point.table),
                        point.column,
                        fmt_plain(point.t9),
                        fmt_float(point.rate),
                        fmt_float(netgen_rate),
                        fmt_float(ratio),
                        f"{log_ratio:.6f}",
                    ]
                )

        stored_count = 0
        stored_bad = 0
        stored_worst: tuple[float, float, float, float, float] | None = None
        for t9, netgen_rate in zip(block.t9, block.rates):
            article_rate = log_interpolate_pairs(finite_points, t9)
            if article_rate is None or article_rate <= 0.0 or netgen_rate <= 0.0:
                continue
            stored_count += 1
            log_ratio = math.log10(netgen_rate / article_rate)
            ratio = netgen_rate / article_rate
            candidate = (abs(log_ratio), t9, article_rate, netgen_rate, log_ratio, ratio)
            if stored_worst is None or candidate[0] > stored_worst[0]:
                stored_worst = candidate
            if abs(log_ratio) > log_threshold:
                stored_bad += 1
                stored_mismatches.append(
                    [
                        reaction,
                        str(first_point.table),
                        first_point.column,
                        fmt_plain(t9),
                        fmt_float(article_rate),
                        fmt_float(netgen_rate),
                        fmt_float(ratio),
                        f"{log_ratio:.6f}",
                    ]
                )

        summary_rows.append(
            [
                reaction,
                str(first_point.table),
                first_point.column,
                "matched",
                str(n_article_finite),
                str(runtime_count),
                str(runtime_bad),
                "..." if runtime_worst is None else f"{runtime_worst[0]:.6f}",
                "..." if runtime_worst is None else fmt_plain(runtime_worst[1]),
                "..." if runtime_worst is None else fmt_float(runtime_worst[2]),
                "..." if runtime_worst is None else fmt_float(runtime_worst[3]),
                "..." if runtime_worst is None else fmt_float(runtime_worst[5]),
                str(stored_count),
                str(stored_bad),
                "..." if stored_worst is None else f"{stored_worst[0]:.6f}",
                "..." if stored_worst is None else fmt_plain(stored_worst[1]),
                "..." if stored_worst is None else fmt_float(stored_worst[2]),
                "..." if stored_worst is None else fmt_float(stored_worst[3]),
                "..." if stored_worst is None else fmt_float(stored_worst[5]),
            ]
        )

    summary_path = output_dir / "iliadis2001_netgen_comparison_summary.dat"
    with summary_path.open("w") as handle:
        handle.write("# Source tables: ppn_nova/references/iliadis2001_table03.dat through table09.dat\n")
        handle.write("# Netgen table: netgen_iliadis2001_log_100.txt\n")
        handle.write("# runtime_* compares PPN-style netgen interpolation at article-listed T9 values.\n")
        handle.write("# stored_* compares stored netgen grid values to log-log interpolation of article rates.\n")
        handle.write(f"# Mismatch threshold: factor > {factor_threshold:g}\n")
        handle.write(
            "\t".join(
                [
                    "reaction",
                    "table",
                    "article_column",
                    "status",
                    "article_finite_points",
                    "runtime_compared_points",
                    "runtime_factor_threshold_mismatches",
                    "runtime_max_abs_log10_ratio",
                    "runtime_worst_T9_GK",
                    "runtime_worst_article_rate",
                    "runtime_worst_netgen_rate",
                    "runtime_worst_netgen_over_article",
                    "stored_compared_points",
                    "stored_factor_threshold_mismatches",
                    "stored_max_abs_log10_ratio",
                    "stored_worst_T9_GK",
                    "stored_worst_article_interp_rate",
                    "stored_worst_netgen_rate",
                    "stored_worst_netgen_over_article",
                ]
            )
            + "\n"
        )
        for row in summary_rows:
            handle.write("\t".join(row) + "\n")

    mismatch_header = [
        "reaction",
        "table",
        "article_column",
        "T9_GK",
        "article_rate",
        "netgen_rate",
        "netgen_over_article",
        "log10_netgen_over_article",
    ]
    runtime_path = output_dir / "iliadis2001_netgen_runtime_mismatches.dat"
    with runtime_path.open("w") as handle:
        handle.write("# PPN-style netgen interpolation at article-listed T9 values.\n")
        handle.write(f"# Rows shown differ by factor > {factor_threshold:g}.\n")
        handle.write("\t".join(mismatch_header) + "\n")
        for row in sorted(runtime_mismatches, key=lambda item: abs(float(item[-1])), reverse=True):
            handle.write("\t".join(row) + "\n")

    stored_path = output_dir / "iliadis2001_netgen_stored_grid_mismatches.dat"
    with stored_path.open("w") as handle:
        handle.write("# Stored netgen grid values compared to article log-log interpolation.\n")
        handle.write(f"# Rows shown differ by factor > {factor_threshold:g}.\n")
        handle.write("\t".join(mismatch_header) + "\n")
        for row in sorted(stored_mismatches, key=lambda item: abs(float(item[-1])), reverse=True):
            handle.write("\t".join(row) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, default=default_pdf_path())
    parser.add_argument("--out", type=Path, default=default_output_dir())
    parser.add_argument("--netgen", type=Path, default=default_netgen_path())
    parser.add_argument("--factor-threshold", type=float, default=2.0)
    args = parser.parse_args()

    lines = extract_pdf_layout(args.pdf)
    tables = parse_tables(lines)
    write_article_tables(tables, args.out, args.pdf)

    netgen = parse_netgen(args.netgen)
    write_comparisons(article_series(tables), netgen, args.out, args.factor_threshold)

    print(f"Wrote Iliadis 2001 Tables 3-9 and netgen comparison files to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

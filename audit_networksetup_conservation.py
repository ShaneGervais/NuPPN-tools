#!/usr/bin/env python3
"""Audit a NuPPN networksetup.txt file for reaction topology problems.

The audit separates hard A/Z rule failures from softer reaction-type
signature mismatches. A signature mismatch means the printed projectile or
ejectile does not agree with PPN's integer reaction type.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable


RTYPE_LABELS = {
    1: "(n,g)",
    2: "(g,n)",
    3: "(n,p)",
    4: "(n,a)",
    5: "(p,g)",
    6: "(g,p)",
    7: "(p,n)",
    8: "(p,a)",
    9: "(a,g)",
    10: "(g,a)",
    11: "(a,n)",
    12: "(a,p)",
    13: "(-,g)",
    14: "(+,g)",
    15: "(b,n)",
    16: "(b,p)",
    17: "(b,a)",
    18: "(v,-)",
    19: "(v,+)",
    20: "(v,p)",
    21: "(v,n)",
    22: "(v,a)",
}

EXPECTED_SIDE_PARTICLES = {
    1: ((1, "NEUT"), (0, "OOOOO")),
    2: ((0, "OOOOO"), (1, "NEUT")),
    3: ((1, "NEUT"), (1, "PROT")),
    4: ((1, "NEUT"), (1, "HE  4")),
    5: ((1, "PROT"), (0, "OOOOO")),
    6: ((0, "OOOOO"), (1, "PROT")),
    7: ((1, "PROT"), (1, "NEUT")),
    8: ((1, "PROT"), (1, "HE  4")),
    9: ((1, "HE  4"), (0, "OOOOO")),
    10: ((0, "OOOOO"), (1, "HE  4")),
    11: ((1, "HE  4"), (1, "NEUT")),
    12: ((1, "HE  4"), (1, "PROT")),
    13: ((0, "OOOOO"), (0, "OOOOO")),
    14: ((0, "OOOOO"), (0, "OOOOO")),
    15: ((0, "OOOOO"), (1, "NEUT")),
    16: ((0, "OOOOO"), (1, "PROT")),
    17: ((0, "OOOOO"), (1, "HE  4")),
    18: ((0, "OOOOO"), (0, "OOOOO")),
    19: ((0, "OOOOO"), (0, "OOOOO")),
    20: ((0, "OOOOO"), (1, "PROT")),
    21: ((0, "OOOOO"), (1, "NEUT")),
    22: ((0, "OOOOO"), (1, "HE  4")),
}

EXPECTED_TOTAL_DELTAS = {
    **{i: (0, 0) for i in range(1, 13)},
    13: (0, 1),
    14: (0, -1),
    15: (0, 1),
    16: (0, -1),
    17: (0, 0),
    18: (0, 1),
    19: (0, -1),
    20: (0, 0),
    21: (0, 0),
    22: (0, 0),
}

SPECIES_RE = re.compile(
    r"^\s*(?P<idx>\d+)\s+(?P<name>.{5})\s+(?P<active>[TF])\s+"
    r"(?P<a>[+-]?\d+(?:\.\d*)?)\s+(?P<z>[+-]?\d+(?:\.\d*)?)\s+"
    r"(?P<state>\d+)\s*$"
)

REACTION_RE = re.compile(
    r"^\s*(?P<idx>\d+)\s+(?P<active>[TF])\s+"
    r"(?P<c1>\d+)\s+(?P<s1>.{5})\s+\+\s+"
    r"(?P<c2>\d+)\s+(?P<s2>.{5})\s+->\s+"
    r"(?P<c3>\d+)\s+(?P<s3>.{5})\s+\+\s+"
    r"(?P<c4>\d+)\s+(?P<s4>.{5})\s+"
    r"(?P<rate>\S+)\s+(?P<source>\S+)\s+(?P<rtype>\S+)\s+"
    r"(?P<ilabb>\d+)\s+(?P<rfac>\S+)\s+(?P<q>\S+)\s*$"
)


@dataclass(frozen=True)
class Species:
    name: str
    a: int
    z: int


@dataclass(frozen=True)
class Reaction:
    line_no: int
    index: int
    active: str
    coeffs: tuple[int, int, int, int]
    species: tuple[str, str, str, str]
    rate: str
    source: str
    rtype: str
    ilabb: int
    rfac: str
    q: str
    text: str


@dataclass(frozen=True)
class Issue:
    line_no: int
    index: int
    active: str
    severity: str
    category: str
    source: str
    rtype: str
    ilabb: int
    lhs: str
    rhs: str
    lhs_a: int
    lhs_z: int
    rhs_a: int
    rhs_z: int
    delta_a: int
    delta_z: int
    expected: str
    note: str


def norm_species(raw: str) -> str:
    return raw.rstrip()


def compact_species(name: str) -> str:
    stripped = name.strip()
    if stripped in {"OOOOO", "PROT", "NEUT"}:
        return stripped
    return re.sub(r"\s+", "", stripped)


def term(coeff: int, species: str) -> str:
    if coeff == 0:
        return "gamma"
    label = compact_species(species)
    return label if coeff == 1 else f"{coeff}{label}"


def reaction_sides(row: Reaction) -> tuple[str, str]:
    c1, c2, c3, c4 = row.coeffs
    s1, s2, s3, s4 = row.species
    return f"{term(c1, s1)} + {term(c2, s2)}", f"{term(c3, s3)} + {term(c4, s4)}"


def parse_networksetup(path: Path) -> tuple[dict[str, Species], list[Reaction]]:
    species: dict[str, Species] = {
        "OOOOO": Species("OOOOO", 0, 0),
    }
    reactions: list[Reaction] = []

    for line_no, line in enumerate(path.read_text().splitlines(), start=1):
        species_match = SPECIES_RE.match(line)
        if species_match:
            name = norm_species(species_match.group("name"))
            species[name] = Species(
                name=name,
                a=round(float(species_match.group("a"))),
                z=round(float(species_match.group("z"))),
            )
            continue

        reaction_match = REACTION_RE.match(line)
        if reaction_match:
            coeffs = tuple(int(reaction_match.group(f"c{i}")) for i in range(1, 5))
            names = tuple(norm_species(reaction_match.group(f"s{i}")) for i in range(1, 5))
            reactions.append(
                Reaction(
                    line_no=line_no,
                    index=int(reaction_match.group("idx")),
                    active=reaction_match.group("active"),
                    coeffs=coeffs,  # type: ignore[arg-type]
                    species=names,  # type: ignore[arg-type]
                    rate=reaction_match.group("rate"),
                    source=reaction_match.group("source"),
                    rtype=reaction_match.group("rtype"),
                    ilabb=int(reaction_match.group("ilabb")),
                    rfac=reaction_match.group("rfac"),
                    q=reaction_match.group("q"),
                    text=line,
                )
            )

    return species, reactions


def species_totals(
    coeffs: Iterable[int], names: Iterable[str], species: dict[str, Species]
) -> tuple[int, int, list[str]]:
    a_total = 0
    z_total = 0
    missing: list[str] = []
    for coeff, name in zip(coeffs, names):
        data = species.get(name)
        if data is None:
            missing.append(name)
            continue
        a_total += coeff * data.a
        z_total += coeff * data.z
    return a_total, z_total, missing


def side_particle_matches(row: Reaction, expected: tuple[tuple[int, str], tuple[int, str]]) -> bool:
    (expected_c2, expected_s2), (expected_c4, expected_s4) = expected
    c1, c2, c3, c4 = row.coeffs
    s1, s2, s3, s4 = row.species

    projectile_ok = c2 == expected_c2 and s2 == expected_s2
    if not projectile_ok and expected_c2 > 0:
        projectile_ok = c2 == 0 and s2 == "OOOOO" and s1 == expected_s2 and c1 > expected_c2

    ejectile_ok = c4 == expected_c4 and s4 == expected_s4
    if not ejectile_ok and expected_c4 > 0:
        ejectile_ok = c4 == 0 and s4 == "OOOOO" and s3 == expected_s4 and c3 > expected_c4

    return projectile_ok and ejectile_ok


def audit_reactions(species: dict[str, Species], reactions: list[Reaction]) -> list[Issue]:
    issues: list[Issue] = []

    for row in reactions:
        lhs, rhs = reaction_sides(row)
        c1, c2, c3, c4 = row.coeffs
        s1, s2, s3, s4 = row.species
        lhs_a, lhs_z, missing_lhs = species_totals((c1, c2), (s1, s2), species)
        rhs_a, rhs_z, missing_rhs = species_totals((c3, c4), (s3, s4), species)
        delta_a = rhs_a - lhs_a
        delta_z = rhs_z - lhs_z
        missing = missing_lhs + missing_rhs

        if missing:
            issues.append(
                Issue(
                    row.line_no,
                    row.index,
                    row.active,
                    "ERROR",
                    "unknown_species",
                    row.source,
                    row.rtype,
                    row.ilabb,
                    lhs,
                    rhs,
                    lhs_a,
                    lhs_z,
                    rhs_a,
                    rhs_z,
                    delta_a,
                    delta_z,
                    "all species present in isotope table",
                    "Missing species: " + ", ".join(compact_species(s) for s in missing),
                )
            )
            continue

        expected_label = RTYPE_LABELS.get(row.ilabb)
        if expected_label is not None and row.rtype != expected_label:
            issues.append(
                Issue(
                    row.line_no,
                    row.index,
                    row.active,
                    "WARN",
                    "rtype_label_mismatch",
                    row.source,
                    row.rtype,
                    row.ilabb,
                    lhs,
                    rhs,
                    lhs_a,
                    lhs_z,
                    rhs_a,
                    rhs_z,
                    delta_a,
                    delta_z,
                    expected_label,
                    "Printed reaction label does not match the integer reaction type.",
                )
            )

        expected_delta = EXPECTED_TOTAL_DELTAS.get(row.ilabb)
        if expected_delta is None:
            if (delta_a, delta_z) != (0, 0):
                issues.append(
                    Issue(
                        row.line_no,
                        row.index,
                        row.active,
                        "ERROR",
                        "custom_type_conservation_violation",
                        row.source,
                        row.rtype,
                        row.ilabb,
                        lhs,
                        rhs,
                        lhs_a,
                        lhs_z,
                        rhs_a,
                        rhs_z,
                        delta_a,
                        delta_z,
                        "Delta A=0, Delta Z=0 for unclassified/custom rows",
                        "Unknown reaction type; checked only total nuclear A and Z.",
                    )
                )
        elif (delta_a, delta_z) != expected_delta:
            category = "weak_rule_violation" if row.ilabb >= 13 else "strong_em_conservation_violation"
            issues.append(
                Issue(
                    row.line_no,
                    row.index,
                    row.active,
                    "ERROR",
                    category,
                    row.source,
                    row.rtype,
                    row.ilabb,
                    lhs,
                    rhs,
                    lhs_a,
                    lhs_z,
                    rhs_a,
                    rhs_z,
                    delta_a,
                    delta_z,
                    f"Delta A={expected_delta[0]}, Delta Z={expected_delta[1]}",
                    "Total nuclear A/Z change does not match the reaction type.",
                )
            )

        expected_sides = EXPECTED_SIDE_PARTICLES.get(row.ilabb)
        if expected_sides is not None and not side_particle_matches(row, expected_sides):
            expected_projectile, expected_ejectile = expected_sides
            issues.append(
                Issue(
                    row.line_no,
                    row.index,
                    row.active,
                    "WARN",
                    "rtype_signature_mismatch",
                    row.source,
                    row.rtype,
                    row.ilabb,
                    lhs,
                    rhs,
                    lhs_a,
                    lhs_z,
                    rhs_a,
                    rhs_z,
                    delta_a,
                    delta_z,
                    f"projectile={term(*expected_projectile)}, ejectile={term(*expected_ejectile)}",
                    "Printed projectile/ejectile does not match PPN's integer reaction type.",
                )
            )

        if row.ilabb in {1, 5, 9} and c3 != 1:
            issues.append(
                Issue(
                    row.line_no,
                    row.index,
                    row.active,
                    "WARN",
                    "capture_label_multiparticle_product",
                    row.source,
                    row.rtype,
                    row.ilabb,
                    lhs,
                    rhs,
                    lhs_a,
                    lhs_z,
                    rhs_a,
                    rhs_z,
                    delta_a,
                    delta_z,
                    "capture rows should have one residual nucleus plus gamma",
                    "The row is labelled as a capture but the residual product coefficient is not 1.",
                )
            )

    return issues


def write_tsv(path: Path, issues: list[Issue]) -> None:
    header = [
        "line_no",
        "index",
        "active",
        "severity",
        "category",
        "source",
        "rtype",
        "ilabb",
        "lhs",
        "rhs",
        "lhs_A",
        "lhs_Z",
        "rhs_A",
        "rhs_Z",
        "delta_A",
        "delta_Z",
        "expected",
        "note",
    ]
    rows = ["\t".join(header)]
    for issue in issues:
        rows.append(
            "\t".join(
                str(value)
                for value in (
                    issue.line_no,
                    issue.index,
                    issue.active,
                    issue.severity,
                    issue.category,
                    issue.source,
                    issue.rtype,
                    issue.ilabb,
                    issue.lhs,
                    issue.rhs,
                    issue.lhs_a,
                    issue.lhs_z,
                    issue.rhs_a,
                    issue.rhs_z,
                    issue.delta_a,
                    issue.delta_z,
                    issue.expected,
                    issue.note,
                )
            )
        )
    path.write_text("\n".join(rows) + "\n")


def markdown_issue_table(issues: list[Issue], limit: int) -> str:
    if not issues:
        return "No rows in this category.\n"
    lines = [
        "| line | index | active | source | type | reaction | delta | note |",
        "|---:|---:|:---:|:---|:---|:---|:---|:---|",
    ]
    for issue in issues[:limit]:
        reaction = f"{issue.lhs} -> {issue.rhs}"
        delta = f"dA={issue.delta_a}, dZ={issue.delta_z}"
        note = issue.note.replace("|", "/")
        lines.append(
            f"| {issue.line_no} | {issue.index} | {issue.active} | {issue.source} | "
            f"{issue.rtype}/{issue.ilabb} | {reaction} | {delta} | {note} |"
        )
    if len(issues) > limit:
        lines.append(f"\nShowing {limit} of {len(issues)} rows. See the TSV for the full list.")
    return "\n".join(lines) + "\n"


def write_report(
    path: Path,
    network_path: Path,
    tsv_path: Path,
    species: dict[str, Species],
    reactions: list[Reaction],
    issues: list[Issue],
    table_limit: int,
) -> None:
    issue_counts = Counter(issue.category for issue in issues)
    active_issue_counts = Counter(issue.category for issue in issues if issue.active == "T")
    severity_counts = Counter(issue.severity for issue in issues)
    active_reactions = sum(1 for row in reactions if row.active == "T")

    hard_categories = {
        "unknown_species",
        "strong_em_conservation_violation",
        "weak_rule_violation",
        "custom_type_conservation_violation",
    }
    hard_issues = [issue for issue in issues if issue.category in hard_categories]
    active_hard = [issue for issue in hard_issues if issue.active == "T"]
    active_signature = [
        issue
        for issue in issues
        if issue.active == "T"
        and issue.category
        in {"rtype_label_mismatch", "rtype_signature_mismatch", "capture_label_multiparticle_product"}
    ]

    ne22_rows = [
        row
        for row in reactions
        if any(name.strip() == "NE 22" for name in row.species)
        and any(name.strip() == "PROT" for name in row.species)
    ]
    ne22_lines = []
    for row in ne22_rows[:8]:
        lhs, rhs = reaction_sides(row)
        ne22_lines.append(
            f"- line {row.line_no}, index {row.index}, active={row.active}: "
            f"{lhs} -> {rhs}, source={row.source}, type={row.rtype}/{row.ilabb}"
        )
    if not ne22_lines:
        ne22_lines.append("- No NE22/proton rows found.")

    summary_lines = [
        "| category | all rows | active rows |",
        "|:---|---:|---:|",
    ]
    for category in sorted(issue_counts):
        summary_lines.append(
            f"| {category} | {issue_counts[category]} | {active_issue_counts[category]} |"
        )

    active_source_counts = Counter(issue.source for issue in issues if issue.active == "T")
    source_lines = [
        "| source | active issue rows |",
        "|:---|---:|",
    ]
    for source, count in active_source_counts.most_common():
        source_lines.append(f"| {source} | {count} |")

    report = f"""# networksetup.txt Conservation Audit

Input: `{network_path}`

Full issue table: `{tsv_path}`

## Summary

- Parsed isotopes: {len(species)}
- Parsed reactions: {len(reactions)}
- Active reactions: {active_reactions}
- Total issues: {len(issues)}
- Literal A/Z or weak-rule issues: {len(hard_issues)}
- Active literal A/Z or weak-rule issues: {len(active_hard)}
- Warning-level labelling/signature issues: {severity_counts["WARN"]}

{chr(10).join(summary_lines)}

Active issue rows by source:

{chr(10).join(source_lines)}

## Important NE22 Note

`22Ne + p -> 23Na + gamma` is balanced: A changes from 23 to 23 and Z from 11 to 11. The nearby non-literal example in this file is the active `22Ne + n -> 22Ne + p` row, which changes total Z from 10 to 11 as printed.

{chr(10).join(ne22_lines)}

## Active Literal Conservation/Rule Issues

{markdown_issue_table(active_hard, table_limit)}

## Active Reaction-Type Signature Warnings

{markdown_issue_table(active_signature, table_limit)}

## How PPN Uses These Rows

There are two relevant PPN paths.

1. Fixed-network mode (`ININET = 1`) reads the topology directly from `networksetup.txt`. In `physics/source/networksetup.F90`, `read_networksetup` assigns the printed input species to `k1/k3` and printed product species to `k7/k5`. In that mode, an active bad product in `networksetup.txt` becomes the actual solver topology.

2. Generated-network mode (`ININET = 0`) calls `networkI`, then `rnetw2007` builds each reaction from target `k1`, branch product `k7`, and the integer reaction type `ilabb`. The projectile/ejectile are hard-coded from `ilabb` in `physics/source/ppn_physics.F90`.

The rate lookup itself is mostly not keyed by the printed product isotope. `physics/source/rates.F90::rates_hash_locations_for_merge` maps common sources by target A/Z and integer reaction type: Reaclib/JINA uses `reaclib_ntrans(ant,znt,1,ilabb)`, NACRE uses `netgen_nacre_ntrans(ant,znt,1,ilabb)`, and ILI01 uses `netgen_illi_ntrans(ant,znt,1,ilabb)`. So a mislabeled product can still receive the rate for the target/type/source, but the abundance flow will follow the topology in `k7/k5`.

## Network-Boundary Interpretation

Most active literal A/Z failures are generated-network boundary mappings, not isolated typos. `physics/source/network_boundaries.F90::natashamcloane` says reaction targets may be adjusted by up to four instantaneous decays after captures or photodisintegrations. `networkI` calls that routine when assigning Reaclib, NACRE, ILI01, and Starlib product branches.

That explains rows such as `21Na + p -> 22Na + gamma`: the ILI01 source reaction is a proton capture to `22Mg`, but if `22Mg` is outside this network the branch can be folded to the in-network beta-decay daughter `22Na`. The printed row omits the positron/neutrino, so it does not conserve charge as a standalone nuclear equation, but the topology is the one PPN evolves.

This also explains rows like `22Ne + n -> 22Ne + p`: the direct `(n,p)` residual would be `22F`, and the generated topology can fold that residual back to `22Ne` after beta decay. Again, the line is not a literal one-step reaction, but the flow target is intentional under the boundary algorithm.

## Rule Notes

- Strong/electromagnetic particle reactions `(n,g)` through `(a,p)` are checked for total nuclear `Delta A = 0`, `Delta Z = 0`.
- Beta-minus and neutrino electron-capture-like rows are checked for `Delta Z = +1`; beta-plus/electron-capture-like rows for `Delta Z = -1`.
- Beta-delayed neutron rows are checked as beta-minus delayed neutron emission (`Delta Z = +1` including the emitted neutron). Beta-delayed proton rows are checked as beta-plus delayed proton emission (`Delta Z = -1` including the emitted proton).
- Custom/unclassified rows such as `(v,v)`/`99` are checked only for total A/Z conservation.
- Signature warnings check whether the printed projectile/ejectile agree with PPN's integer reaction type. They allow PPN's collapsed notation where identical reactants or products are combined.
"""
    path.write_text(report)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "networksetup",
        nargs="?",
        default="ppn_nova/networksetup.txt",
        type=Path,
        help="networksetup.txt file to audit",
    )
    parser.add_argument(
        "--tsv",
        type=Path,
        default=Path("ppn_nova/references/networksetup_conservation_issues.tsv"),
        help="output TSV file",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("ppn_nova/references/networksetup_conservation_report.md"),
        help="output Markdown report",
    )
    parser.add_argument(
        "--table-limit",
        type=int,
        default=75,
        help="maximum number of rows per Markdown issue table",
    )
    args = parser.parse_args()

    species, reactions = parse_networksetup(args.networksetup)
    issues = audit_reactions(species, reactions)

    args.tsv.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.tsv, issues)
    write_report(args.report, args.networksetup, args.tsv, species, reactions, issues, args.table_limit)

    hard_categories = {
        "unknown_species",
        "strong_em_conservation_violation",
        "weak_rule_violation",
        "custom_type_conservation_violation",
    }
    hard_issues = [issue for issue in issues if issue.category in hard_categories]
    active_hard = [issue for issue in hard_issues if issue.active == "T"]
    print(f"parsed reactions: {len(reactions)}")
    print(f"issues written: {len(issues)}")
    print(f"literal conservation/rule issues: {len(hard_issues)}")
    print(f"active literal conservation/rule issues: {len(active_hard)}")
    print(f"wrote: {args.tsv}")
    print(f"wrote: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

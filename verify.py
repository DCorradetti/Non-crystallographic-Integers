#!/usr/bin/env python3
"""Top-level verification driver for the nc-integers certificate archive.

Re-runs each of the six protocol steps described in the article's
*Computational methods* section, regenerates the certificates under
``certificates/``, and prints a SHA-256 digest of the certificate tree.

This driver assumes the scripts in ``scripts/`` are self-contained and
write their certificates into the directories under ``certificates/``.

No third-party dependencies. Python >= 3.11.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts"
CERTS = ROOT / "certificates"

STEPS: list[tuple[str, str, list[str], str]] = [
    (
        "P1+P2",
        "Closure, integrality, and H_4 shell of the icosian double",
        ["golden_octonions.py", "--out", str(CERTS / "icosian_double")],
        "icosian_double",
    ),
    (
        "P3+P4",
        "Bounded G3 search and Gram-matrix self-duality",
        ["g3_search.py", "--out", str(CERTS / "g3")],
        "g3",
    ),
    (
        "P3",
        "Dual-discriminant certificate (det G_I = phi^2, det G_{G_0} = 2+3 phi)",
        ["g3_dual_discriminant.py", "--out", str(CERTS / "g3")],
        "g3",
    ),
    (
        "P3 cont.",
        "Candidate denominator-ideal sweep against the dual-discriminant",
        ["g3_candidate_ideals_from_discriminant.py", "--cert-dir", str(CERTS / "g3")],
        "g3",
    ),
    (
        "P5",
        "Ramified sqrt(5)-denominator scan over F_5^8",
        ["g3_denominator_sqrt5_scan.py", "--out", str(CERTS / "g3")],
        "g3",
    ),
    (
        "P3 audit",
        "E_8 + E_8 unimodularity / 5-modular trace lattice audit",
        ["g3_e8e8_unimodularity.py", "--out", str(CERTS / "g3")],
        "g3",
    ),
    (
        "P6",
        "Trace-integral G3-B' search (h4 and sqrt5 ansatz)",
        ["g3_trace_integral_search.py", "--out", str(CERTS / "g3_b_prime")],
        "g3_b_prime",
    ),
    (
        "P6 cont.",
        "G3-B' phase 2: octonion-stable isotropic subspaces in (Z/5Z)^8",
        ["g3_b_prime_phase2.py", "--out", str(CERTS / "g3_b_prime")],
        "g3_b_prime",
    ),
]


def _run(label: str, description: str, cmd: list[str]) -> None:
    script = SCRIPTS / cmd[0]
    if not script.exists():
        print(f"[{label}] ERROR: script not found: {script}", file=sys.stderr)
        sys.exit(2)
    print(f"\n=== [{label}] {description} ===")
    print(f"--> python {script.name} {' '.join(cmd[1:])}")
    proc = subprocess.run(
        [sys.executable, str(script), *cmd[1:]],
        cwd=str(SCRIPTS),
        check=False,
    )
    if proc.returncode != 0:
        print(f"[{label}] FAILED with exit code {proc.returncode}", file=sys.stderr)
        sys.exit(proc.returncode)


def _hash_tree(root: Path) -> str:
    """Deterministic SHA-256 of a directory tree.

    Hashes the sorted relative paths together with the SHA-256 of each file.
    Writes the per-file digests into ``certificates/SHA256SUMS``.
    """
    digest = hashlib.sha256()
    sumlines: list[str] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.name == "SHA256SUMS":
            continue
        rel = path.relative_to(root).as_posix()
        file_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        sumlines.append(f"{file_hash}  {rel}")
        digest.update(rel.encode())
        digest.update(b"\0")
        digest.update(bytes.fromhex(file_hash))
    sums_path = root / "SHA256SUMS"
    sums_path.write_text("\n".join(sumlines) + "\n", encoding="utf-8")
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove certificates/ before regenerating (recommended for verification).",
    )
    parser.add_argument(
        "--only",
        type=str,
        default=None,
        help="Comma-separated list of step labels to run (default: all).",
    )
    args = parser.parse_args()

    if args.clean and CERTS.exists():
        for sub in CERTS.iterdir():
            if sub.is_dir():
                shutil.rmtree(sub)
        for f in CERTS.glob("*"):
            if f.is_file():
                f.unlink()

    CERTS.mkdir(exist_ok=True)
    for sub in ("icosian_double", "g3", "g3_b_prime"):
        (CERTS / sub).mkdir(exist_ok=True)

    only = set(s.strip() for s in args.only.split(",")) if args.only else None
    for label, description, cmd, _subdir in STEPS:
        if only is not None and label not in only:
            print(f"-- skipping [{label}]")
            continue
        _run(label, description, cmd)

    digest = _hash_tree(CERTS)
    print("\n========================================")
    print(f"Certificate tree SHA-256: {digest}")
    print(f"Per-file digests written to: {CERTS / 'SHA256SUMS'}")
    print("========================================")


if __name__ == "__main__":
    main()

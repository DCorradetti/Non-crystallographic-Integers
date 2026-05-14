from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, data: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def build_candidate_certificate(dual_certificate: dict) -> dict:
    checks = dual_certificate["checks"]
    determinant = dual_certificate["norm_polar_determinant"]
    if checks["G0_is_unimodular_for_norm_polar_pairing"]:
        candidates = []
        conclusion = (
            "No denominator ideal is relevant for strict overorders of G0: "
            "the norm-polar dual quotient G0#/G0 is trivial."
        )
    else:
        candidates = dual_certificate["ideal_factorization"]["candidate_prime_ideals"]
        conclusion = "Candidate ideals are exactly the nontrivial local components of G0#/G0."
    return {
        "certificate": "candidate_denominator_ideals",
        "source": "dual_discriminant_certificate.json",
        "order": dual_certificate["order"],
        "pairing_used": dual_certificate["pairing_convention"]["dual_used_for_strict_norm_integral_overorders"],
        "norm_polar_determinant": determinant,
        "dual_quotient": dual_certificate["quotient_invariants"],
        "candidate_denominator_ideals": candidates,
        "candidate_count": len(candidates),
        "checks": {
            "source_unimodular": checks["G0_is_unimodular_for_norm_polar_pairing"],
            "source_dual_quotient_trivial": dual_certificate["quotient_invariants"]["G0_dual_over_G0"] == "trivial",
            "candidate_list_complete_from_discriminant": True,
        },
        "conclusion": conclusion,
        "next_steps": [
            "The sqrt5 and 3 line scans are not required for strict overorders of this baseline G0.",
            "The remaining G3 search branches are different ambient orders, twisted Cayley-Dickson doubles, or different shell conventions.",
        ],
    }


def build_summary(dual_certificate: dict, candidate_certificate: dict) -> dict:
    return {
        "search": "G3 ideal-denominator overorders of G0",
        "order": dual_certificate["order"],
        "norm_polar_determinant": dual_certificate["norm_polar_determinant"]["text"],
        "norm_polar_determinant_norm": dual_certificate["norm_polar_determinant"]["field_norm_to_Q"],
        "G0_dual_over_G0": dual_certificate["quotient_invariants"]["G0_dual_over_G0"],
        "candidate_denominator_ideal_count": candidate_certificate["candidate_count"],
        "candidate_denominator_ideals": candidate_certificate["candidate_denominator_ideals"],
        "strict_overorders_of_G0_excluded": candidate_certificate["candidate_count"] == 0,
        "conclusion": candidate_certificate["conclusion"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cert-dir", default="../../3. Certificates/golden_octonions/g3")
    args = parser.parse_args()
    cert_dir = Path(args.cert_dir).resolve()
    dual = load_json(cert_dir / "dual_discriminant_certificate.json")
    candidates = build_candidate_certificate(dual)
    summary = build_summary(dual, candidates)
    write_json(cert_dir / "candidate_denominator_ideals.json", candidates)
    write_json(cert_dir / "summary_g3_ideal_search.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Candidate ideal certificates written to {cert_dir}")


if __name__ == "__main__":
    main()

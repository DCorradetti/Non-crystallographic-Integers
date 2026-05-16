from __future__ import annotations

import argparse
import json
from pathlib import Path

from nc_quasicrystal_core import CERT_DIR, h2_artifacts, write_json


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(CERT_DIR))
    args = parser.parse_args()
    out = Path(args.out)
    artifacts = h2_artifacts()
    for name, data in artifacts.items():
        write_json(out / name, data)
    print(json.dumps({"written": sorted(artifacts), "out": str(out)}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

from __future__ import annotations

import json

from nc_quasicrystal_core import run_all


def main() -> None:
    run_dir = run_all()
    manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

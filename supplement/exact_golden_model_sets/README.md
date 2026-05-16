# Exact Golden Model Sets — Supplementary Archive

Reproducibility package for

> D. Corradetti, *Exact golden arithmetic of H2 and H4 model sets*.

The package contains the Python scripts that regenerate the certificates
listed in the paper's Reproducibility Manifest, the regenerated JSON
certificates themselves, and the Markdown audit notes.

## Layout

```
supplement/exact_golden_model_sets/
  README.md
  scripts/
    nc_quasicrystal_core.py            (canonical core module)
    golden_octonions.py                (bundled dependency)
    script_exact_golden_model_sets_v2.py
    script_h2_model_set_v2.py
    script_h4_projection_v2.py
    script_k8_e8_reciprocal_v2.py
  certificates/
    Certificate_SG1_h2_benchmark_v2.md
    Certificate_SG2_h4_projection_v1.md
    Certificate_SG3_k8_e8_reciprocal_v2.md
    Certificate_SG4_reproducibility_v1.md
    icosian_trace_gram.json            (bundled dependency)
    summary.json
    run_manifest.json
    cut_project_basis.json
    diffraction_intensities_h2.json
    h2_shell_certificate.json
    h4_shell_coordinates.json
    icosahedral_star_vectors.json
    k8_e8_projection_comparison.json
    penrose_patch_vertices.json
    projection_matrices_h4.json
    reciprocal_module_h2.json
    reciprocal_module_h4.json
```

## Requirements

- Python >= 3.11
- sympy

## Running the Full Verification

From the repository root:

```bash
python supplement/exact_golden_model_sets/scripts/script_exact_golden_model_sets_v2.py
```

This runs the full pipeline and writes regenerated JSON certificates to
`supplement/exact_golden_model_sets/certificates/`, plus a timestamped
immutable copy to `supplement/exact_golden_model_sets/runs/`.

## Focused Wrappers

Each wrapper regenerates one slice of the certificate package:

```bash
python supplement/exact_golden_model_sets/scripts/script_h2_model_set_v2.py
python supplement/exact_golden_model_sets/scripts/script_h4_projection_v2.py
python supplement/exact_golden_model_sets/scripts/script_k8_e8_reciprocal_v2.py
```

## Path Auto-Detection

`nc_quasicrystal_core.py` auto-detects two layouts:

- the author's local working tree (`2. Scripts/`, `3. Certificates/`,
  `7. Results/`), and
- this self-contained supplement layout (`scripts/`, `certificates/`,
  `runs/`).

The bundled `golden_octonions.py` and `icosian_trace_gram.json` make the
supplement runnable without the rest of the working tree.

## Scope

The package supports the mathematical model-set and reciprocal-module
results of the paper. It makes no experimental materials claim and does
not use the octonionic associator data of the companion articles.

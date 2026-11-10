import argparse, hashlib, json, os, platform, sys, time
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
ROOT = Path(__file__).resolve().parent

def pipeline():
    ap = argparse.ArgumentParser(description='Reproduce the grouped-emission computational case study.')
    ap.add_argument('--reps', type=int, default=1000, help='Number of replications per simulation scenario (default: 1000)')
    args = ap.parse_args()
    if args.reps < 2:
        ap.error('--reps must be at least 2')

    required = [ROOT / 'data' / 'grouped_inputs.csv', ROOT / 'data' / 'resampled_assignments.csv']
    required += [ROOT / 'data' / f'reference_table_{i:02}.csv' for i in range(1, 12)]
    missing = [str(p.relative_to(ROOT)) for p in required if not p.is_file()]
    if missing:
        border = "=" * 72
        msg = (
            f"\n{border}\n"
            "Supplementary Input Files Required\n"
            f"{border}\n"
            "The numerical inputs used in this study are presented in Table 1 of the\n"
            "article. Supplementary Excel files (samplecounts.xlsx and 18-cell_data.xlsx)\n"
            "are provided separately through the journal (see DATA_AVAILABILITY.md).\n\n"
            "To execute the computational pipeline, place the required input files in\n"
            "the 'data/' directory. Currently missing:\n"
            + "\n".join(f"  - {m}" for m in missing)
            + f"\n{border}\n\n"
        )
        sys.stderr.write(msg)
        sys.exit(1)

    sys.path.insert(0, str(ROOT / 'src'))
    from analyze import main
    from extension import run
    from verify import verify
    from figures import figures
    from tables import tables
    from report import report

    start = time.perf_counter()
    main(str(ROOT / 'data' / 'grouped_inputs.csv'), str(ROOT / 'data' / 'resampled_assignments.csv'), str(ROOT / 'results'), args.reps)
    run(ROOT, args.reps)
    tables(ROOT)
    verify(ROOT, args.reps)
    figures(ROOT)
    report(ROOT)

    import numpy, scipy, openpyxl, matplotlib
    versions = {k: m.__version__ for k, m in [('numpy', numpy), ('scipy', scipy), ('openpyxl', openpyxl), ('matplotlib', matplotlib)]}
    versions.update(
        python=platform.python_version(),
        platform=platform.platform(),
        replications=args.reps,
        elapsed_seconds=round(time.perf_counter() - start, 3)
    )
    (ROOT / 'validation' / 'execution.json').write_text(json.dumps(versions, indent=2))
    print('Pipeline completed. See validation/verification.json and results/manuscript_comparison.csv.', flush=True)

    paths = [p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and '.git' not in p.parts and p.name not in ['SHA256_MANIFEST.json', '.DS_Store']]
    (ROOT / 'SHA256_MANIFEST.json').write_text(json.dumps({str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}, indent=2))

if __name__ == '__main__':
    pipeline()

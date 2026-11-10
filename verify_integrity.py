import hashlib, json
from pathlib import Path

root = Path(__file__).resolve().parent
manifest_path = root / 'SHA256_MANIFEST.json'
if not manifest_path.is_file():
    raise SystemExit('Error: SHA256_MANIFEST.json not found.')

manifest = json.loads(manifest_path.read_text())
failures = []
for name, expected in manifest.items():
    path = root / name
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        failures.append(name)

if failures:
    raise SystemExit('Integrity check failed: ' + ', '.join(failures))

print(f'Integrity verified for {len(manifest)} files.')

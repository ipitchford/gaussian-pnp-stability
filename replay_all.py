#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Replay every exact verifier from a fresh temporary copy of this package and
compare the regenerated records with the shipped ones.

Needs sympy==1.13.3 (with mpmath) for certify_revision.py, emit_certificates.py
and practical_numbers.py; the other scripts use only the standard library.
Fields that record timing or interpreter versions are excluded from comparison.
The figure (make_figure.py) is illustration only and is not replayed here.
Exit status 0 means every script passed and every compared record matched.
"""
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = ['verify_exact.py', 'verify_circulant.py', 'verify_reentrant.py', 'certify_revision.py',
           'emit_certificates.py', 'check_certificates.py', 'practical_numbers.py', 'make_fingerprints.py']
SAME_BYTES = ['certificates/thm4_classification.json', 'certificates/thm4_tstar.json',
              'certificates/thm3_circulant.json', 'circulant_certificate.json', 'exact_fingerprints.json']
VOLATILE = {'elapsed_seconds', 'python', 'sympy'}


def strip(obj):
    if isinstance(obj, dict):
        return {k: strip(v) for k, v in obj.items() if k not in VOLATILE}
    if isinstance(obj, list):
        return [strip(v) for v in obj]
    return obj


def main():
    bad = []
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp) / 'pkg'
        shutil.copytree(HERE, work, ignore=shutil.ignore_patterns('.git', '__pycache__'))
        for name in SCRIPTS:
            r = subprocess.run([sys.executable, '-I', '-B', name], cwd=work, capture_output=True, text=True)
            print(f'{"PASS" if r.returncode == 0 else "FAIL"}: {name} (exit {r.returncode})')
            if r.returncode:
                bad.append(name)
                print(r.stdout[-2000:], r.stderr[-2000:])
            if name == 'check_certificates.py' and r.stdout != (HERE / 'certificate_check_output.txt').read_text():
                bad.append('check_certificates.py stdout differs from certificate_check_output.txt')
        for rel in SAME_BYTES:
            same = (work / rel).read_bytes() == (HERE / rel).read_bytes()
            print(f'{"PASS" if same else "FAIL"}: {rel} byte-identical')
            if not same:
                bad.append(rel)
        for rel in ['revision_certificate.json', 'reentrant_verification.json']:
            same = strip(json.loads((work / rel).read_text())) == strip(json.loads((HERE / rel).read_text()))
            print(f'{"PASS" if same else "FAIL"}: {rel} identical apart from timing/version fields')
            if not same:
                bad.append(rel)
    print('REPLAY PASSED' if not bad else f'REPLAY FAILED: {bad}')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()

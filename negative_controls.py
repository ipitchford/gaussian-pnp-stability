#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Semantic negative controls for the Appendix C certificate checker.

Each control copies certificates/ to a temporary directory, applies one
plausible corruption and requires check_certificates.py to FAIL. The unmodified
copy must PASS first. Standard library only. Exit status 0 means every
corruption was rejected.
"""
import json
import shutil
import subprocess
import sys
import tempfile
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHECK = HERE / 'check_certificates.py'


def load(d, n):
    return json.loads((d / n).read_text())


def save(d, n, obj):
    (d / n).write_text(json.dumps(obj, indent=1) + '\n')


def c_bernstein_sign(d):
    c = load(d, 'thm4_classification.json')
    e = c['side_conditions'][0]
    e['numerator_bernstein']['coefficients'][0][0] = str(-Fraction(e['numerator_bernstein']['coefficients'][0][0]))
    save(d, 'thm4_classification.json', c)


def c_side_condition_poly(d):
    c = load(d, 'thm4_classification.json')
    e = c['side_conditions'][1]
    i, j, v = e['numerator'][0]
    e['numerator'][0] = [i, j, str(-Fraction(v) * 10**6)]
    save(d, 'thm4_classification.json', c)


def c_delta_coefficient(d):
    c = load(d, 'thm4_classification.json')
    i, j, v = c['Delta'][0]
    c['Delta'][0] = [i, j, str(Fraction(v) * 2)]
    save(d, 'thm4_classification.json', c)


def c_eps_interval(d):
    c = load(d, 'thm4_classification.json')
    c['eps_crit_interval'] = ['1/50000', '1/40000']
    save(d, 'thm4_classification.json', c)


def c_tstar_interval(d):
    c = load(d, 'thm4_tstar.json')
    c['t_star_interval'] = ['53/100', '531/1000']
    save(d, 'thm4_tstar.json', c)


def c_tstar_coefficient(d):
    c = load(d, 'thm4_tstar.json')
    c['integer_coeffs_ascending'][1] = str(int(c['integer_coeffs_ascending'][1]) // 2)
    save(d, 'thm4_tstar.json', c)


def c_tc_interval(d):
    c = load(d, 'thm3_circulant.json')
    c['t_c_interval'] = ['19157074942/100000000000000', '19157074943/100000000000000']
    save(d, 'thm3_circulant.json', c)


def c_tc_constant_term(d):
    c = load(d, 'thm3_circulant.json')
    c['integer_coeffs_ascending'][0] = str(-int(c['integer_coeffs_ascending'][0]))
    save(d, 'thm3_circulant.json', c)


CONTROLS = [c_bernstein_sign, c_side_condition_poly, c_delta_coefficient, c_eps_interval,
            c_tstar_interval, c_tstar_coefficient, c_tc_interval, c_tc_constant_term]


def run(d):
    return subprocess.run([sys.executable, '-I', '-B', str(CHECK), '--dir', str(d)],
                          capture_output=True, text=True).returncode


def main():
    bad = []
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp) / 'baseline'
        shutil.copytree(HERE / 'certificates', base)
        ok = run(base) == 0
        print(f'{"PASS" if ok else "FAIL"}: unmodified certificates accepted')
        if not ok:
            bad.append('baseline')
        for k, f in enumerate(CONTROLS):
            d = Path(tmp) / f'c{k}'
            shutil.copytree(HERE / 'certificates', d)
            f(d)
            rejected = run(d) != 0
            print(f'{"PASS" if rejected else "FAIL"}: {f.__name__[2:]} rejected')
            if not rejected:
                bad.append(f.__name__)
    print('ALL NEGATIVE CONTROLS PASSED' if not bad else f'NEGATIVE CONTROLS FAILED: {bad}')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Reproduce the reviewed C/N/S marker panel using verified external inputs."""
from pathlib import Path
import argparse
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
import time

SECTION = Path(__file__).resolve().parent
REPO = SECTION.parents[1]
sys.path.insert(0, str(REPO/'checks'))
from check_functional_cycles_inputs import digest, verify


def run_step(script, output, env, arguments=()):
    started = time.monotonic()
    command = [sys.executable, '-u', str(script), *map(str, arguments)]
    print('START '+script.name, flush=True)
    with (output/'logs'/(script.stem+'.log')).open('w') as log:
        child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=env)
        try:
            while True:
                try:
                    code = child.wait(timeout=30)
                    break
                except subprocess.TimeoutExpired:
                    print(f'RUNNING {script.name}: {time.monotonic()-started:.0f} s', flush=True)
        except BaseException:
            child.terminate()
            try:
                child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait()
            raise
    if code:
        print((output/'logs'/(script.stem+'.log')).read_text()[-6000:], file=sys.stderr)
        raise RuntimeError(f'{script.name} failed with exit code {code}')
    print(f'PASS {script.name}: {time.monotonic()-started:.1f} s', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--with-previous-comparison', action='store_true',
                        help='Include the descriptive prior-analysis and denominator audit; needs three extra inputs.')
    parser.add_argument('--calibrate', action='store_true',
                        help='Also repeat design simulations, including the rejected bootstrap-t interval diagnostic.')
    parser.add_argument('--no-figures', action='store_true')
    parser.add_argument('--reference-results', type=Path,
                        help='Optionally compare generated numeric tables with an external reviewed result directory.')
    args = parser.parse_args()
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        parser.error('Do not disable Python assertions; they guard the scientific producer scripts.')
    output = args.output_dir.resolve()
    if output.exists():
        parser.error('Output directory already exists; choose a fresh directory.')
    if args.reference_results and not args.reference_results.is_dir():
        parser.error('Reference-results directory does not exist.')
    versions = {p: importlib.metadata.version(p) for p in ['numpy', 'pandas', 'scipy', 'matplotlib']}
    try:
        preflight, entries = verify(args.input_dir, args.with_previous_comparison)
    except ValueError as error:
        parser.error(str(error))
    output.mkdir(parents=True, exist_ok=False)
    for name in ['inputs', 'results', 'figures', 'logs']:
        (output/name).mkdir()
    shutil.copytree(SECTION/'code', output/'code', ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(SECTION/'analysis_specification.json', output/'analysis_specification.json')
    for item in entries:
        destination = output/'inputs'/item['path']
        shutil.copy2(args.input_dir/item['path'], destination)
        if digest(destination) != item['sha256']:
            raise RuntimeError('Input changed while staging: '+item['path'])
    manifest = {'status': 'RUNNING', 'input_check': preflight, 'input_files': entries,
                'source_code': {p.name: digest(p) for p in sorted((output/'code').glob('*.py'))},
                'specification_sha256': digest(output/'analysis_specification.json'),
                'python': platform.python_version(), 'packages': versions,
                'calibration_requested': args.calibrate, 'figures_requested': not args.no_figures,
                'reference_comparison_requested': bool(args.reference_results)}
    record = output/'execution_manifest.json'
    record.write_text(json.dumps(manifest, indent=2)+'\n')
    env = dict(os.environ)
    for name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        env[name] = '1'
    env['MPLCONFIGDIR'] = str(output/'matplotlib_cache')
    steps = []
    if args.calibrate:
        steps += ['calibrar_y_verificar.py', 'verificar_intervalos_HC3.py']
    steps += ['ejecutar_cns.py', 'verificar_resultados.py', 'subunidades_apoyo.py', 'report_tables.py']
    if args.with_previous_comparison:
        steps += ['comparar_analisis_previos.py']
    if not args.no_figures:
        steps += ['figuras_cns.py']
    try:
        for name in steps:
            run_step(output/'code'/name, output, env)
        validate_args = ['--run-root', output]
        if args.reference_results:
            validate_args += ['--reference-results', args.reference_results.resolve()]
        run_step(REPO/'checks/verify_functional_cycles.py', output, env, validate_args)
        manifest['status'] = 'PASS'
    except BaseException:
        manifest['status'] = 'FAILED'
        record.write_text(json.dumps(manifest, indent=2)+'\n')
        raise
    manifest['output_files'] = {str(p.relative_to(output)): digest(p)
                                for directory in ['results', 'figures']
                                for p in sorted((output/directory).iterdir()) if p.is_file()}
    record.write_text(json.dumps(manifest, indent=2)+'\n')
    print('COMPLETE: '+str(output), flush=True)


if __name__ == '__main__':
    main()

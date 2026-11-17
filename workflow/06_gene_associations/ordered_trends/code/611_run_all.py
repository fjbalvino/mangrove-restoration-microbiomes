#!/usr/bin/env python3
"""Execute the complete analysis in Python; R replica is independent."""
from pathlib import Path
import argparse, importlib.util, traceback

def load(name):
    spec=importlib.util.spec_from_file_location(name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();out=Path(a.out)
    if out.exists():raise SystemExit('Output directory already exists. Use a new directory.')
    try:
        load('611_monotonic').run(a.inputs,a.out)
        (out/'STATUS.txt').write_text('RUNNING_SENSITIVITY_LOLO\n')
        load('611_sensitivity_LOLO').run(a.inputs,a.out)
        (out/'STATUS.txt').write_text('RUNNING_FIGURES\n')
        load('611_figures').run(a.inputs,a.out)
        (out/'STATUS.txt').write_text('COMPLETE\n')
    except Exception:
        if out.exists():(out/'STATUS.txt').write_text('FAILED\n'+traceback.format_exc())
        raise

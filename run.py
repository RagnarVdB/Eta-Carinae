#!/usr/bin/env python

# Gebruik:
# ./run.py par_bestand.par mod_usr.t output_map

import argparse
import re
from os import getcwd, path
from pathlib import Path
from shutil import rmtree
from subprocess import run

THREADS = 4
WORKDIR = Path(path.dirname(path.realpath(__file__))) / "workdir"

if not WORKDIR.exists():
    WORKDIR.mkdir()

current_dir = Path(getcwd())

parser = argparse.ArgumentParser(description="Run MPI-AMRVAC")
parser.add_argument("par_file", type=str, help="parameter file")
parser.add_argument("mod_usr", type=str, help="mod_usr file")
parser.add_argument("output_dir", type=str, help="map voor .vtu files")
parser.add_argument(
    "--purge",
    action=argparse.BooleanOptionalAction,
    help="verwijder bestaande bestanden in output_dir",
)
args = parser.parse_args()


out_dir = args.output_dir

x = Path(out_dir)

if x.exists() and args.purge:
    rmtree(x)

if not x.exists():
    x.mkdir()


with open(args.par_file, "r") as f:
    text = f.read()
    regex = r"(?<=base_filename=')(.*)(?=')"
    old = re.findall(regex, text)[0]
    if "/" in old:
        folder, base = old.split("/")
        new = str(current_dir / out_dir / base)
    else:
        new = str(current_dir / out_dir / old)
    new_file = re.sub(regex, new, text)

with open(WORKDIR / "par_file.par", "w") as f:
    f.write(new_file)


run(f"cp {args.mod_usr} {WORKDIR/'mod_usr.t'}", shell=True)
run(f"cd {WORKDIR} && setup.pl -d=2", shell=True)
run(f"cd {WORKDIR} && make clean && make", shell=True)
run(f"cd {WORKDIR} && mpirun -np {THREADS} ./amrvac -i par_file.par", shell=True)

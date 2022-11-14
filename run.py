#!/usr/bin/env python

# Gebruik:
# ./run.py mijn_par_bestand.par mijn_mod_usr.t mijn_output_map


# workdir instellen
# mkdir workdir
# verander WORKDIR in python file

# Script executable:
# chmod 744 run.py


import re
from os import getcwd, path
from pathlib import Path
from subprocess import run
from sys import argv

THREADS = 4
WORKDIR = Path(path.dirname(path.realpath(__file__))) / "workdir"
# Map waarin AMRVAC runt

current_dir = Path(getcwd())

arguments = argv[1:]

if len(arguments) == 3:
    par_file, mod_usr, out_dir = arguments

    with open(par_file, "r") as f:
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

else:
    raise ValueError("2 of 3 argumenten verwacht")


run(f"cp {mod_usr} {WORKDIR/'mod_usr.t'}", shell=True)
run(f"cd {WORKDIR} && setup.pl -d=2", shell=True)
run(f"cd {WORKDIR} && make clean && make", shell=True)
run(f"cd {WORKDIR} && mpirun -np {THREADS} ./amrvac -i par_file.par", shell=True)

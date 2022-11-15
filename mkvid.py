#!/usr/bin/env python

# Gebruik:
# ./mkvid par_bestand.par output_map videonaam
# of
# ./mkvid par_bestand.par output_map -b base_bestandsnaam -f framerate videonaam

import argparse
from python_imaging import create_video
from pathlib import Path
from os import getcwd

current = Path(getcwd())

parser = argparse.ArgumentParser(description="Maak video van amrvac output")
parser.add_argument("par_file", type=str, help="parameter file")
parser.add_argument("output_folder", type=str, help="map met .vtu files")
parser.add_argument("-v", "--variable", type=str, help="Geplotte variabele 'rho', 'p', 'v1', 'v2' (rho default)")
parser.add_argument("-f", "--framerate", type=int, help="framerate")
parser.add_argument("-b", "--base", type=str, help="basenaam vtu bestanden")
parser.add_argument("output_video", help="bestandsnaam output video", type=str)

args = parser.parse_args()

par_file = Path(args.par_file)
assert par_file.exists()
output_folder = Path(args.output_folder)
assert output_folder.exists()
output_video = Path(args.output_video)
if args.variable is None:
    variable = "rho"
else:
    variable = args.variable

base = args.base

if args.framerate is None:
    framerate = 1
else:
    framerate = args.framerate

if base is None:
    iss = []
    base = None
    for file in output_folder.glob("*.vtu"):
        new_base = file.stem[:-4]
        number = file.stem[-4:]
        if new_base != base and base is not None:
            raise ValueError("More than one base in folder")
        base = new_base
        i = int(number)
        iss.append(i)
else:
    iss = []
    for file in output_folder.glob(f"{base}*.vtu"):
        i = int(file.name[len(base) : -4])
        iss.append(i)

create_video(
    sorted(iss),
    str(par_file),
    base,
    str(output_folder) + "/",
    variable,
    framerate,
    str(current/output_video),
)

# Project Eta-Cariae

## Paper:
https://www.aanda.org/articles/aa/pdf/2022/03/aa42020-21.pdf

## Varnames
https://amrvac.org/md_doc_varnames.html

## Scriptjes:
executable maken:

    chmod 744 run.py && chmod 744 mkvid.py

### Help:
    ./run.py --help
    ./mkvid.py --help
### Running
    ./run.py par_bestand.par mod_usr_bestand.t output_map --purge

of
```
setup.pl -d=2
make clean
make -j 4
mpirun -np 4 ./amrvac -i amrvac.par
```

### Video maken
    ./mkvid par_bestand.par output_map videonaam
of

    ./mkvid par_bestand.par output_map -b base_bestandsnaam -f framerate video_naam
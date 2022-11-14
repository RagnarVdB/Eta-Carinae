# Project Eta-Cariae

## Paper:
https://www.aanda.org/articles/aa/pdf/2022/03/aa42020-21.pdf

## Varnames
https://amrvac.org/md_doc_varnames.html

## Running
```bash
./run.py par_bestand.par mod_usr_bestand.t output_map
```

```bash
setup.pl -d=2
make clean
make -j 4
mpirun -np 4 ./amrvac -i amrvac.par
```

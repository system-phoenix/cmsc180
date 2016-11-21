#!/bin/bash

mpicc -o lab06.exe lab06.c -lm
sbatch submit.sh ./lab06.exe
cat output.txt

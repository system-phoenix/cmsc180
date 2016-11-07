#!/bin/bash

mpicc -o sample.exe sample.c
sbatch submit.sh ./sample.exe
cat output.txt

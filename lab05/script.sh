#!/bin/bash

mpicc -o lab05.exe lab05.c
sbatch submit.sh ./lab05.exe
cat output.txt

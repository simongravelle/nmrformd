#!/bin/bash
#OAR -n nmrformd
#OAR -l /nodes=1/cpu=1/core=1,walltime=48:00:00
#OAR --stdout log.out
#OAR --stderr log.err
#OAR --project tamtam

source /bettik/gravells/venv/bin/activate

python3 run_nmrformd.py


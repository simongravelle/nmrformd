#!/bin/bash

cd references/
python3 filter-reference.py
cd ..

make clean
make html

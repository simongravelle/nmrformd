#!/bin/bash

git submodule update --remote

cd references/
python3 filter-reference.py
cd ..

make clean
make html

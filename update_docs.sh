#!/bin/bash
shopt -s extglob

cd docs
rm -rf -- !("build"|"source")
cd ../

make clean
make html

shopt -s dotglob
mv docs/build/html/* docs/
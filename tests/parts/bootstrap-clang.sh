#!/bin/bash

rm -rf build

CXX=clang++ CC=clang cmake -Bbuild -H. -DWITH_TCMALLOC=On -DCMAKE_EXPORT_COMPILE_COMMANDS=1 

mv build/compile_commands.json .

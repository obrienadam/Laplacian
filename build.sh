#!/usr/bin/env bash

mpic++ laplacian.cc -o laplacian -lmpi -lpetsc -std=c++11 -O3 -ffast-math -Wl,-rpath=/usr/local/lib

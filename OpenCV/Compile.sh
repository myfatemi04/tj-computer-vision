#!/bin/bash
g++ $1.cpp -o $1 -Wall `pkg-config --cflags --libs opencv4`

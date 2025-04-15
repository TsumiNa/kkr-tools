#!/bin/bash

#====================================================
# sym_prog: program path to execute calculation.
# 
# This variable must be setted.
#====================================================
sym_prog="../sym-program/Release/bin/sym"


#======================================================
# execute script
#=====================================================
ln -sf ${sym_prog} sym

./sym > @Al-sym2

#!/bin/bash

rm -f @Al-sym2
rm -f symm_coef
rm -f sym
rm -f harmony_input
rm -f lattice_input
rm -f harmony_output
rm -f gf_input

if [ -n "$1" ] && [ "$1" == "all" ]; then
    rm -f *.bak
    rm -f *.coeff
    rm -f shape*
    rm -f greenout_win
fi

#!/bin/bash

#====================================================
# tbkkr_prog: program path to execute calculation.
#
# This variable must be setted.
#====================================================
tbkkr_prog="../tbkkr-program/Release/bin/bulk"


#=================================================
# OMP_NUM_THREADS, MKL_NUM_THREADS
# define number of thread used by parallel calculation
# if you don't know what it is, keep default
#================================================
export MKL_DYNAMIC=true
# export OMP_NUM_THREADS=4
# export MKL_NUM_THREADS=8
unset OMP_NUM_THREADS
unset MKL_NUM_THREADS

# lsda calculation
lsda () {
    cp inputcard-lsda inputcard
    ./tbkkr > @Pd-lsda
    rm -f inputcard

    cp potin.740 potin.740.save &> /dev/null
    cp fort.23  gf_indata &> /dev/null
    cp fort.11  potin.740 &> /dev/null
    cp fort.3   gen.pot &> /dev/null

    rm -f  fort.* &> /dev/null
}

# gga calculation
gga () {
    cp  inputcard-gga  inputcard
    ./tbkkr > @Pd-gga
    rm -f  inputcard

 
    cp  fort.23 gf_indata &> /dev/null
    cp  fort.11 potin.740.green &> /dev/null
    cp  fort.3 gen.pot &> /dev/null
    rm -f  fort.* &> /dev/null
}
 

#======================================================
# execute script
#=====================================================
ln -sf ${tbkkr_prog} tbkkr
rm -f \ green &> /dev/null

if [ -n "$1" ]; then
    if [ "$1" == "lsda" ]; then
        lsda
    elif [ "$1" == "gga" ]; then
        gga
    else
        echo "paramter must be [lsda] or [gga]."
        echo "none paramter for run lsda then gga sequentially."
    fi
else
    lsda
    gga
fi

rm -f \ green
rm -f tbkkr


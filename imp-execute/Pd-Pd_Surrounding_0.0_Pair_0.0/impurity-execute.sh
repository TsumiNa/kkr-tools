#!/bin/bash

#====================================================
# bulk_root: set the path of green functions and tu.coeff, back.coeff file.
# green: which green you want to use for this calculation.
# imp_prog: program path to execute calculation.
# 
# These variables must be setted
#====================================================
bulk_root="../../sym-execute"
green="greenout_win"
imp_prog="../../imp-program/Release/bin/impurity"


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


#======================================================
# init
#=====================================================
init () {
    rm -f green &> /dev/null
    rm -f *.coeff &> /dev/null
    rm -f imp &> /dev/null

    ln -sf ${bulk_root}/${green} green
    ln -sf ${bulk_root}/back.coeff back.coeff
    ln -sf ${bulk_root}/tu.coeff tu.coeff
    ln -sf ${imp_prog} imp    
}


#======================================================
# cleanup
#=====================================================
cleanup () {
    rm -f fort.* &> /dev/null
    rm -f wfct* &> /dev/null
    rm -f *.BAK &> /dev/null
    rm -f broy &> /dev/null
    rm -f green_TRANS &> /dev/null
}


#======================================================
# execute lsda
#=====================================================
lsda () {
    ./imp < input_relax-lsda > @Pd102l4-lsda
    mv fort.11 potout-lsda &> /dev/null
    for pot in potential-* do
    do
        mv ${pot} ${pot}.save
        cp potout-lsda-scf ${pot}
    done

    cleanup
}



#======================================================
# execute gga
#=====================================================
gga () {
    ./imp < input_relax-nsgga > @Pd102l4-nsgga
    cp fort.11 potout-nsgga
    cleanup
}


#======================================================
# execute script
#=====================================================
if [ -n "$1" ]; then
    if [ "$1" == "lsda" ]; then
        init
        lsda
    elif [ "$1" == "gga" ]; then
        init
        gga
    else
        echo "paramter must be [lsda] or [gga]."
        echo "none paramter for run lsda then gga sequentially."
    fi
else
    init
    lsda
    gga
fi

rm -f green &> /dev/null
rm -f *.coeff &> /dev/null
rm -f imp &> /dev/null

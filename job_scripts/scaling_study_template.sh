#!/bin/bash

#Place any necessary job scheduler commands above

##############################################################################
# System-specific commands (e.g., module load, copying rayleigh etc.) go here


##############################################################################
# Rayleigh scaling section
mpicmd="mpiexec -np"  # The mpiexec or mpirun command up to the process count
raexec="rayleigh.opt" # The rayleigh executable to test

maxcpu=2048 # Won't test process combinations with ncpu > maxcpu
mincpu=1024 # and similarly so for ncpu < mincpu

nmn=2 
nmx=7 # Will test values of nprow and npcol from 2^nmn up to 2^nmx

# Problem size specification
# Nr must be even.  Suggest that nr is a power of 2
# Ntheta should be 3 x a power of 2.
# each pair of nrs[i], nts[i] will be tested
# The setup below defines four problems sizes (not 16)

nrs=( "64"  "128" "256" "512" )  # nr values 
nts=( "192" "384" "768" "1536" ) # ntheta values

# If desired, 
istart=1
iend=1



###############################################################
# You should not need to edit below this line

i=$istart
while [ $i -le $iend ] # Begin outer loop over problem sizes
do

    nr=${nrs[$i]}
    colmax=$nr

    nt=${nts[$i]}
    rowmax=$((nt / 3))

    n=$nmn
    while [ $n -le $nmx ] # Begin loop over nprow values
    do
        nprow=$((2 ** $n))
      
        m=$nmn
        while [ $m -le $nmx ]  # Begin loop over npcol values
        do

            npcol=$((2 ** $m))
            ncpu=$(($nprow * $npcol))
            
            if [ $npcol -le $colmax ] && [ $nprow -le $rowmax ] && [ $ncpu -le $maxcpu ] && [ $ncpu -ge $mincpu ]
            then

            runcmd=$mpicmd" "$ncpu" "$raexec" -nprow "$nprow" -npcol "$npcol" -nr "$nr" -nt "$nt
            echo $runcmd
            
            fi
            
            m=$(($m + 1))
        done                  # End loop over npcol values
            
        n=$(($n + 1))

    done                      # End loop over nprow values

    i=$(($i + 1))

done    # End Loop over problem sizes



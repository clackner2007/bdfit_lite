#!/usr/bin/env bash
#################
# Finds the condor outputs that did not work by checking which condor/out* files did not write output files
#
##################

#for condor
#grep 'writing output file' condor/out* | sed 's/condor\/out\([^:]*\):.*/shell\1.sh/' > done
#ls shell*.sh > attempt #| sed 's/bigRun13113\///' > attempt
#rm -f d2
#for f in `cat attempt`; do grep $f done >> d2; done
#diff d2 attempt | grep '^>' | sed 's/> //' > failed
#rm -f done attempt
#echo 'failed shell scripts listed in failed'
#rm -f attempt done

nm=`grep "#PBS -N" $1 | awk '{print $3}'`
echo $nm
grep 'writing output file' pbs/*$nm*out* | sed 's/.*out-\([0-9]\+\):.*/\1/' > done
grep "#PBS -t" $1 | grep -o [0-9]* > attempt
comm -13 <(sort done) <(sort attempt) | sort -n > failed

#!/usr/bin/env python
#############
# make_condor_shell.py
#
#############
"""
creates the shells needed to submit condor jobs to the queue
To actually submit:
"""

import argparse
import os


def writeShell(infile, filenum, numperjob, initgal=0, outdir='./', imagedir='./',
               profiles=r'{DVC:8}'):
    """
    writes a shell file to wrap each IDL job
    """
    shellfile = open('shell{0:04d}.sh'.format(filenum), 'w')
    runline = "echo \"fit_sample, '{0}', {1:d}, {2:d}, '{3}', '{4}', profiles={5}, /residuals, /cutoff\" | idl\n".\
    format(infile, filenum * numperjob + initgal, (filenum+1) * numperjob + initgal,
           outdir, imagedir, profiles)
    shellfile.write(r"""#!/usr/bin/env bash
###
# shell script to run fit_sample.pro
####
export SHELL=/bin/bash
source /u/dss/products/eups/bin/setups.sh
setup idlutils
setup idl #get the correct version
source ${HOME}/.idlenv

""")
    shellfile.write(runline+"\n")
    shellfile.close()
    
    os.chmod('shell{0:04d}.sh'.format(filenum), 0755)
    

def writeSubmit(filenum):
    """
    writes a submit file for each shell file
    """
    subfile = open('submit.{0:04d}'.format(filenum), 'w')
    subfile.write(r"""universe=vanilla
requirements=(Arch =="X86_64" && OpSys == "LINUX")
executable=shell{0:04d}.sh
error=condor/err{0:04d}.sh
log=/tmp/condor/log{0:04d}.sh
output=condor/out{0:04d}
queue 1
""".format(filenum))
    subfile.close()


def writeDAG(njobs, maxjobs=15):
    """
    writes the condor DAG file which allows multiple job submissions at once
    (like a scatter-gather)
    """
    dagfile = open('dag_submit', 'w')
    
    for j in range(njobs):
        dagfile.write('JOB job{0:04d} submit.{0:04d}\n'.format(j))
    for j in range(njobs):
        dagfile.write('CATEGORY job{0:04d} groups\n'.format(j))
    dagfile.write('MAXJOBS groups {0:d}\n'.format(maxjobs))
    dagfile.close()    
    

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile', help='input .FITS file')
    parser.add_argument('numgal', type=int, help='total number of galaxies')
    parser.add_argument('numperJob', type=int, help='number of galaxies per job')
    parser.add_argument('-f', '--first', dest='initgal', default=0, help='first galaxy in list',
                        type=int)
    parser.add_argument('-o', '--outdir', dest='outdir', default='./', help='output directory')
    parser.add_argument('-i', '--imagedir', dest='imagedir', default='./', help='image directory')
    parser.add_argument('-m', '--maxjobs', dest='maxjobs', default=15, type=int,
                        help='max # of condor jobs running at once')
    args = parser.parse_args()
    
    numjobs, leftover = divmod(args.numgal, args.numperJob)
    if leftover > 0:
        numjobs += 1
    
    for i in range(numjobs):
        writeShell(args.infile, i, args.numperJob, args.initgal, args.outdir, args.imagedir)
        writeSubmit(i)
        
    writeDAG(numjobs, args.maxjobs)
    
    

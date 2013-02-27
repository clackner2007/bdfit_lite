#!/usr/bin/env python

#####make the dag script
#######################
"""
%prog num_jobs namefor_stages submit#(0,or 1)
"""

import sys,os,re,glob
import subprocess
import optparse



def main():

    parser = optparse.OptionParser(usage=__doc__)
    parser.add_option("-m","--max_jobs", dest="maxJOB",type=int,
		      default=15,help="max IDL jobs at once")
    
    opts,args = parser.parse_args()

    

    if len(args) != 3:
	parser.print_help()
	sys.exit(1)
	
    njobs = int(args[0])
    filename = args[1]
    submitnum = args[2]


    file_all = open('/tmp/condor/submits/'+filename+'dag_all','w')
    file_all.write("#multiple dag jobs\n")
    
    submits = glob.glob('/u/clackner/work/bulges/bulgeFit/submit*0'+
			submitnum+'[0-9][0-9]*')
    
    #write a dag file for each run
    for n in range(njobs):
	file_dag = open('/tmp/condor/submits/'+filename+'dag_{0:02d}'.format(n), 'w')
	file_dag.write("#single run through "+filename+"\n")
	for i, s in enumerate(submits):
	    file_dag.write("JOB j{0:04d} ".format(i)+s+"\n")
	for i, s in enumerate(submits):
	    file_dag.write("CATEGORY j{0:04d} groups".format(i)+"\n")
	for i, s in enumerate(submits):
	    file_dag.write("RETRY j{0:04d} 2\n".format(i))
	file_dag.write("MAXJOBS groups {0:d}\n".format(opts.maxJOB))
   
	file_dag.close()

    dags = glob.glob("/tmp/condor/submits/"+filename+"dag_[0-9]*")

    for i, d in enumerate(dags):
    
	file_all.write("SUBDAG EXTERNAL  j{0:02d} ".format(i)+d+"\n")
        #print "JOB j{0:02d} ".format(i)+d+".condor.sub"
	file_all.write("SCRIPT POST j{0:02d} /u/clackner/work/bulges/bulgeFit/".format(i)+
		       filename+"process.sh\n")

    for i,d in enumerate(dags[:-1]):
	file_all.write("PARENT j{0:02d} CHILD j{1:02d}\n".format(i,i+1))

    file_all.close()

if __name__=='__main__':
    main()

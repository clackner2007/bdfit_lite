#!/usr/bin/env perl

############
#
# writes shell scripts for submitting sersic fitting to condor
###########

use warnings;
use strict;
use POSIX;

my ($numGal, $numperJob, $initGal) = @ARGV;
die "Usage: $0 totalnumbergalaxies, number per job first galaxy number\n"
    unless $numGal and $numperJob;

my $numjobs = $numGal/$numperJob + 1;


my $filenum = 0;

foreach my $i (0..$numjobs-1)
{
    my $filenum = $i;
    my $file = sprintf( 'shell%04i.sh', $filenum );
	
    open( OUT, ">$file" );
	
    print OUT "#!/usr/bin/env bash\n".
	"#############\n".
	"# Claire Lackner\n# Mar 1 2010\n#\n".
	"# shell script for sersicfitting fake galaxies in idl\n".
	"#  to be used to submit to condor\n".
	"######\n\n";
	
    print OUT "export SHELL=/bin/bash\n";
    print OUT "source /u/dss/products/eups/bin/setups.sh\n";
    print OUT "#setup idlutils\n";
    print OUT "#source /usr/peyton/common/licensed/idl71/bin/idl_setup.bash\n";
    print OUT "#source /u/clackner/.idlenv\n";
#	print OUT "export IDL_PATH=/u/clackner/work/bulges/bulgeBranch/:\$IDL_PATH\n";
	
    my $inputfile='/u/clackner/work/sandbox/fastfakegals/standards.fits';
    printf OUT "echo \"bulge_fakesample, \'%s\', %u, %u\" | idl\n",
	$inputfile, $i*$numperJob+$initGal, ($i+1)*$numperJob+$initGal;
	
    close( OUT );
    
    system("chmod 0755 $file");



#write individual submit files
    $file = sprintf( 'submit.%04i', $filenum );
    open( SUB, ">$file" );

    print SUB "universe=vanilla\n".
	"requirements = (Arch ==\"X86_64\" && OpSys == \"LINUX\")\n";
	
    printf SUB "executable=shell%04i.sh\n", $filenum;
    print SUB "getenv=True\n";
    printf SUB "error=condor/err%04i\n", $filenum;
    printf SUB "log=/tmp/condor/log%04i\n", $filenum;
    printf SUB "output=condor/out%04i\n", $filenum;
    print SUB "queue 1\n\n";
	
    close( SUB );
	
    $filenum++;
	
}


#write DAG submit file
open( DAG, ">dag_submit" );

print DAG "#\n#DAG submit file for jobs that need to be run\n#\n\n";

foreach my $k (0..$numjobs-1)
{
    printf DAG "JOB j%04i submit.%04i\n", $k, $k;
}
foreach my $j (0..$numjobs-1)
{
    printf DAG "CATEGORY j%04i groups\n", $j;
}

print DAG "MAXJOBS groups 8\n";

close(DAG);

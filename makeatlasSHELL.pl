#!/usr/bin/env perl

############
#
# writes shell scripts for submitting sersic fitting to condor
###########

use warnings;
use strict;
use POSIX;

my ($numGal, $numperJob, $initGal, $inputfile, $outputfile, $filenum) = @ARGV;
die "Usage: $0 totalnumbergalaxies, number per job first galaxy number inputfile, output folder starting file number\n"
    unless $numGal and $numperJob and $inputfile and $outputfile;



my $numjobs = $numGal/$numperJob + 1;


#my $filenum = 0;
my $set = (split(/[\/.]/,$inputfile))[-2];
my $file_first = $filenum;
foreach my $i (0..$numjobs-1)
{
    #my $filenum = $i;
    my $file = sprintf( 'shell%04i%s.sh', $filenum,$set );
	
    open( OUT, ">$file" );
	
    print OUT "#!/usr/bin/env bash\n".
	"#############\n".
	"# Claire Lackner\n# Mar 1 2010\n#\n".
	"# shell script for sersicfitting COSMOS galaxies in idl\n".
	"#  to be used to submit to condor\n".
	"######\n\n";
	
    print OUT "export SHELL=/bin/bash\n";
    print OUT "source /u/dss/products/eups/bin/setups.sh\n";
    print OUT "setup idlutils\n";
    print OUT "setup photoop v1_9_11\n";
    print OUT "setup idl v7_1\n";
    print OUT "#source /usr/peyton/common/licensed/idl71/bin/idl_setup.bash\n";
    print OUT "source /u/clackner/.idlenv\n";
    print OUT "export PHOTO_REDUX=/u/dss/redux\n";
    print OUT "export PHOTO_DATA=/u/dss/data\n";
    print OUT "export SPECTRO_DATA=/u/dss/spectro\n";
    #print OUT "echo \$PATH\n";
    #print OUT "echo \$IDL_PATH\n";
    
#	print OUT "export IDL_PATH=/u/clackner/work/bulges/bulgeBranch/:\$IDL_PATH\n";
	
#    my $inputfile='/u/clackner/work/bulges/data/atlassample3/S0111-01594.fits';
#    my $outputfile = '/u/clackner/work/bulges/outputZeroBound/S0111-01594/';
#    my $outputfile = '/u/clackner/work/bulges/reina/output/';
#    my $inputfile='/u/clackner/work/bulges/reina/reina_rc2.fits';
    
#    printf OUT "echo \"bulge_atlassample, \'%s\', %u, %u, \'%s\',do_filters=[\'r\']\" | idl\n",
#    $inputfile, $i*$numperJob+$initGal, ($i+1)*$numperJob+$initGal, $outputfile;
    
#    printf OUT "echo  \"bulge_fakesample, \'%s\', \'%s\', %u, %u\"| idl\n",
#    $inputfile, $outputfile, $i*$numperJob+$initGal, ($i+1)*$numperJob+$initGal,;
    printf OUT "echo \"bulge_atlassample, \'%s\', %u, %u, \'%s\', \'%s\'\" | idl\n",
    $inputfile, $i*$numperJob+$initGal, ($i+1)*$numperJob+$initGal, $outputfile, '/peyton/scr/chimera0/clackner/RM_COSMOS/';
    close( OUT );
    
    system("chmod 0755 $file");



#write individual submit files
    $file = sprintf( 'submit.%04i%s', $filenum,$set );
    open( SUB, ">$file" );

    print SUB "universe=vanilla\n".
	"requirements = (Arch ==\"X86_64\" && OpSys == \"LINUX\")\n";
	
    printf SUB "executable=shell%04i%s.sh\n", $filenum,$set;
#    print SUB "getenv=True\n";
    printf SUB "error=condor/err%04i%s\n", $filenum,$set;
    printf SUB "log=/tmp/condor/log%04i%s\n", $filenum,$set;
    printf SUB "output=condor/out%04i%s\n", $filenum,$set;
    print SUB "queue 1\n\n";
	
    close( SUB );
	
    $filenum++;
	
}


my @lines;
#read in old DAG submit file
if(open( DAG, "<dag_submit" )){
    my $line = <DAG>;
    while( !eof && $line !~ /^CAT/)
    {
	if( $line =~/^JOB/)
	{
	    push @lines, $line;
	}
	$line = <DAG>;
    }
    close(DAG);
}
my $nlines = scalar @lines + $numjobs;
my $kfirst = scalar @lines;
open( DAG, ">dag_submit100");
    
print DAG "#\n#DAG submit file\n#\n";

foreach my $l (0..$kfirst-1)
{
    printf DAG $lines[$l];
}
foreach my $k (0..$numjobs-1)
{
    printf DAG "JOB j%04i submit.%04i%s\n", $k+$kfirst, $k+$file_first, $set;
}
foreach my $j (0..$nlines-1)
{
    printf DAG "CATEGORY j%04i groups\n", $j;
}
print DAG "MAXJOBS groups 15\n";

close(DAG);

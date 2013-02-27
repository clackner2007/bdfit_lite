#!/usr/bin/env perl

############
#
# writes shell scripts for submitting sersic fitting to condor
###########

use warnings;
use strict;
use POSIX;

my ($lowra, $highra, $lowdec, $highdec) = @ARGV;
die "Usage: $0 lowra, highra, lowdec, highdec\n" 
    unless $lowra and $highdec and $highra and $lowdec;

#set the declination increments to 10 minutes
my $num_decs = ceil($highdec) - floor($lowdec);
my $dec_incr = 0.1667;


#die "invalid declination range" unless $num_decs > 0;

my $filenum = 0;

foreach my $i (0..$num_decs-1)
{
    my $dec1 = floor($lowdec) + $i;

    foreach my $j (0..5)
    {
	my $file = sprintf( 'shell%02i.sh', $filenum );
	
	open( OUT, ">$file" );
	
	print OUT "#!/usr/bin/env bash\n".
	    "#############\n".
	    "# Claire Lackner\n# Sept 17 2009\n#\n".
	    "# shell script for sersicfit in idl\n".
	    "#  to be used to submit to condor\n".
	    "######\n\n";
	
	print OUT "export SHELL=/bin/bash\n";
	print OUT "source /u/dss/products/eups/bin/setups.sh\n";
	print OUT "setup idlutils\n";
	print OUT "source /u/clackner/.idlenv\n";
#	print OUT "export IDL_PATH=/u/clackner/work/bulges/bulgeBranch/:\$IDL_PATH\n";
	
	
	my $dec2 = $dec1 + $dec_incr;
    
	printf OUT "echo \"bulgefit, %+.2f, %+.2f, %+.2f, %+.2f, [1.5, 1.0, 0.5, 0.2]\" | idl\n",
	$lowra, $highra, $dec1, $dec2;
	
	$dec1 = $dec2;
	
	close( OUT );
	
	system("chmod 0755 $file");



#write individual submit files
	$file = sprintf( 'submit.%02i', $filenum );
	open( SUB, ">$file" );

	print SUB "universe=vanilla\n".
	    "requirements = (Arch ==\"X86_64\" && OpSys == \"LINUX\")\n";
	
	printf SUB "executable=shell%02i.sh\n", $filenum;
	printf SUB "error=condor/err%02i\n", $filenum;
	printf SUB "log=/tmp/condor/log%02i\n", $filenum;
	printf SUB "output=condor/out%02i\n", $filenum;
	print SUB "queue 1\n\n";
	
	close( SUB );
	
	$filenum++;
	
    }

}
#write DAG submit file
open( DAG, ">dag_submit" );

print DAG "#\n#DAG submit file for jobs that need to be rerun\n#\n\n";

foreach my $i (0..$filenum-1)
{
    printf DAG "JOB j%02i submit.%02i\n", $i, $i;
}
foreach my $i (0..$filenum-1)
{
    printf DAG "CATEGORY j%02i doubles\n", $i;
}

print DAG "MAXJOBS doubles 10\n";

close(DAG);
# now, write the condor submit file

#open( COND, ">submit_shells" );

#foreach my $i (0..$filenum-1)
#{
#    print COND "universe=vanilla\n".
#	"requirements = (Arch ==\"X86_64\" && OpSys == \"LINUX\")\n";
#    
#    printf COND "executable=shell%02i.sh\n", $i;
#    printf COND "error=condor/err%02i\n", $i;
#    printf COND "log=condor/log%02i\n", $i;
#    printf COND "output=condor/out%02i\n", $i;
#    print COND "queue 1\n\n";
#}
#close( COND );

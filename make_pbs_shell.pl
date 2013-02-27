#!/usr/bin/env perl

############
#
# writes shell scripts for submitting sersic fitting to condor
###########

use warnings;
use strict;
use POSIX;

my ($numGal, $numperJob, $initGal, $inputfile, $outputfile, $fileext, $imagepath) = @ARGV;
die "Usage: $0 totalnumbergalaxies, number per job first galaxy number inputfile, output folder, pbs file id, /path/to/images\n"
    unless $numGal and $numperJob and $inputfile and $outputfile and $imagepath;

my $lastGal = $initGal+$numGal;
my $set = (split(/[\/.]/,$inputfile))[-2];
my $job_list = "$initGal";
my $i=$initGal+$numperJob;
while ($i < $lastGal)
{
    $job_list .= ",$i";
    $i += $numperJob;
}
##write the pbs array file
my $pbs_file=sprintf('pbs_array_%s_%s', $set, $fileext);
open (SUB, ">$pbs_file" );

print SUB "#!/usr/bin/env bash\n";
printf SUB "#PBS -t $job_list\n";
printf SUB "#PBS -N %s_%s\n", $set, $fileext;
print SUB "#PBS -o pbs/\${PBS_JOBNAME}.out\n".
    "#PBS -e pbs/\${PBS_JOBNAME}.err\n".
    "#PBS -v DISPLAY=0,PATH=/usr/bin:/bin:/usr/local/bin\n".
     "cd \$PBS_O_WORKDIR\n";
print SUB ". /home/clackner/.idlenv\n";
print SUB "env\n\n";
print SUB "type idl\n";
printf SUB "echo \"bulge_atlassample, \'%s\', \${PBS_ARRAYID}, \$((\${PBS_ARRAYID}+$numperJob)), \'%s\', \'%s\'\"\n", $inputfile, $outputfile, $imagepath;
printf SUB "echo \"bulge_atlassample, \'%s\', \${PBS_ARRAYID}, \$((\${PBS_ARRAYID}+$numperJob)), \'%s\', \'%s\'\" | idl\n", $inputfile, $outputfile, $imagepath;

close( SUB );

system("chmod u+x $pbs_file");

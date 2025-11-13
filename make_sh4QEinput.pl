=b
make qe input files for all strucutres in labelled folders.
You need to use this script in the dir with all dpgen collections (in all_cfgs folder)
perl ../tool_scripts/cfg2QEinput.pl 
=cut
use warnings;
use strict;
use JSON::PP;
use Data::Dumper;
use List::Util qw(min max);
use Cwd;
use POSIX;
use Parallel::ForkManager;
use List::Util qw/shuffle/;


my $currentPath = getcwd();# dir for all scripts
chdir("..");
my $mainPath = getcwd();# main path of Perl4dpgen dir
chdir("$currentPath");

########## source folder you need to assign
#my $sou_dir = "$currentPath/thermo_label/C_AlPNT-T0260-lmpT50to1250-P0-R500000/labelled";#source dir with labelled subfolders
my $sou_dir = "$currentPath/thermo_label/*/labelled";#source dir with labelled subfolders


my $forkNo = 1;#although we don't have so many cores, only for submitting jobs into slurm
my $pm = Parallel::ForkManager->new("$forkNo");
#for your sh file
my %sbatch_para = (
            nodes => 1,#how many nodes for your qe job
            threads => 1,#'$(nproc)',#1,#modify it to 2, 4 if oom problem appears            
            partition => "All",#which partition you want to use
            runPath => "/opt/thermoPW-7-2_intel/bin/pw.x", #qe executable          
            );
#my $onlyCheckcfgNumber = "no";#only check how many labelled cfg files, then decide how 
#many you want to use for DFT,
## maxmium number for dft input files to create, set a super large number for selecting all 
#my $maxDFTfiles = 30;

#you may use relabel script to increase more labelled folders if needed
my @oldsh = `find $sou_dir -type f -name "*.sh" -exec readlink -f {} \\;|sort|grep -v "C_AlPNT-T0260-lmpT50to1250-P0-R500000"`;
for my $i (@oldsh){`rm -f $i`;}

my @oldsout = `find $sou_dir -type f -name "*.sout" -exec readlink -f {} \\;|sort|grep -v "C_AlPNT-T0260-lmpT50to1250-P0-R500000"`;
for my $i (@oldsout){`rm -f $i`;}

my @allQEin = `find $sou_dir -type f -name "*.in" -exec readlink -f {} \\;|sort|grep -v "C_AlPNT-T0260-lmpT50to1250-P0-R500000"`;
map { s/^\s+|\s+$//g; } @allQEin;
print "***all QE input Number: ",scalar @allQEin,"\n";
#my @pathOfAllcfgs;
my $jobNo = 0;
for my $i (@allQEin){
    print "$i\n";
    my $basename = `basename $i`;
    my $dirname = `dirname $i`;
    $basename =~ s/\.in//g; 
    chomp ($basename,$dirname);
    `rm -f $dirname/$basename.sh`;
    $jobNo++;
my $here_doc =<<"END_MESSAGE";
#!/bin/sh
#SBATCH --output=$basename.sout
#SBATCH --job-name=$basename
#SBATCH --nodes=$sbatch_para{nodes}
#SBATCH --partition=$sbatch_para{partition}
#SBATCH --reservation=script_test
##SBATCH --ntasks-per-node=12
##SBATCH --exclude=node23
##SBATCH --nodelist=master
hostname
source /opt/intel/oneapi/setvars.sh
rm -rf pwscf*
node=$sbatch_para{nodes}
threads=$sbatch_para{threads}
processors=\$(nproc)
np=\$((\$node*\$processors/\$threads))
export OMP_NUM_THREADS=\$threads
#the following two are for AMD CPU if slurm chooses for you!!
export MKL_DEBUG_CPU_TYPE=5
export MKL_CBWR=AUTO
export LD_LIBRARY_PATH=/opt/mpich-4.0.3/lib:/opt/intel/oneapi/mkl/latest/lib:\$LD_LIBRARY_PATH
export PATH=/opt/mpich-4.0.3/bin:\$PATH

/opt/mpich-4.0.3/bin/mpiexec -np \$np $sbatch_para{runPath} -in $basename.in
rm -rf pwscf*
perl /opt/qe_perl/QEout_analysis.pl
perl /opt/qe_perl/QEout2data.pl

END_MESSAGE
    #chomp $here_doc;
    print "writing $dirname/$basename.sh\n";
    open(FH, "> $dirname/$basename.sh") or die $!;
    print FH $here_doc;
    close(FH);       
}#  


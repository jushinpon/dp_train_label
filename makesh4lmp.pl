=b
make lmp input files for all strucutres in labelled folders.
You need to use this script in the dir with all dpgen collections (in all_cfgs folder)
perl ../tool_scripts/cfg2lmpinput.pl 
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

my $sou_dir = "$currentPath/thermo_label";#source dir
my $submitJobs = "no";
my %sbatch_para = (
            nodes => 1,#how many nodes for your lmp job
           # ntasks_per_node => 24,
            partition => "All",#which partition you want to use
            LAMMPSPath => "/opt/lammps-mpich-4.0.3/lmpdeepmd", #lmp executable          
            mpiPath => "/opt/mpich-4.0.3/bin/mpiexec" #mpipath          
            );

my $forkNo = 1;#although we don't have so many cores, only for submitting jobs into slurm
my $pm = Parallel::ForkManager->new("$forkNo");

my @alllmpin = `find $sou_dir -type f -name "*.in" -exec readlink -f {} \\;|sort`;
map { s/^\s+|\s+$//g; } @alllmpin;

my $jobNo = 0;
for my $i (@alllmpin){
    my $basename = `basename $i`;
    $basename =~ s/^\s+|\s+$//g;
    $basename =~ s/\.in//g; 

    my $dirname = `dirname $i`;
    $dirname =~ s/^\s+|\s+$//g;

    `rm -f $dirname/$basename.sh`;
    $jobNo++;
my $here_doc =<<"END_MESSAGE";
#!/bin/sh
#SBATCH --output=$basename.log
#SBATCH --job-name=label_$basename
#SBATCH --nodes=$sbatch_para{nodes}
#SBATCH --partition=$sbatch_para{partition}
hostname
rm -rf *.cfg
rm -rf *.data
rm -rf lmp_output
node=1
threads=2
processors=\$(nproc)
np=\$((\$node*\$processors/\$threads))
export OMP_NUM_THREADS=\$threads
export TF_INTRA_OP_PARALLELISM_THREADS=\$np
export TF_INTER_OP_PARALLELISM_THREADS=\$threads
#The following are used only for intel MKL (not workable for AMD)
export KMP_AFFINITY=granularity=fine,compact,1,0
export KMP_BLOCKTIME=0
export KMP_SETTINGS=TRUE


export LD_LIBRARY_PATH=/opt/mpich-4.0.3/lib:\$LD_LIBRARY_PATH
export PATH=/opt/mpich-4.0.3/bin:\$PATH
export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64_lin:\$LD_LIBRARY_PATH

$sbatch_para{mpiPath} -np \$np $sbatch_para{LAMMPSPath} -in $basename.in

END_MESSAGE

    open(FH, "> $dirname/$basename.sh") or die $!;
    print FH $here_doc;
    close(FH);
    if($submitJobs eq "yes"){
        chdir($dirname);
        `sbatch $basename.sh`;
        chdir($currentPath);
    }    
}#  


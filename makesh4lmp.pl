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
            LAMMPSPath => "lmp", #lmp executable          
            #LAMMPSPath => "/opt/lammps-mpich-4.0.3/lmpdeepmd", #lmp executable          
            mpiPath => "mpiexec" #mpipath          
            #mpiPath => "/opt/mpich-4.0.3/bin/mpiexec" #mpipath          
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
##SBATCH --reservation=script_test
hostname
rm -rf *.cfg
rm -rf *.data
rm -rf lmp_output

if [ -f /opt/anaconda3/bin/activate ]; then
    
    source /opt/anaconda3/bin/activate deepmd-cpu-v3
    export LD_LIBRARY_PATH=/opt/deepmd-cpu-v3/lib:/opt/deepmd-cpu-v3/lib/deepmd_lmp:\$LD_LIBRARY_PATH
    export PATH=/opt/deepmd-cpu-v3/bin:\$PATH

elif [ -f /opt/miniconda3/bin/activate ]; then
    source /opt/miniconda3/bin/activate deepmd-cpu-v3
    export LD_LIBRARY_PATH=/opt/deepmd-cpu-v3/lib:/opt/deepmd-cpu-v3/lib/deepmd_lmp:\$LD_LIBRARY_PATH
    export PATH=/opt/deepmd-cpu-v3/bin:\$PATH
else
    echo "Error: Neither /opt/anaconda3/bin/activate nor /opt/miniconda3/bin/activate found."
    exit 1  # Exit the script if neither exists
fi

node=1
threads=\$(nproc)
processors=\$(nproc)
np=\$((\$node*\$processors/\$threads))

export OMP_NUM_THREADS=\$processors
export TF_INTRA_OP_PARALLELISM_THREADS=\$processors
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


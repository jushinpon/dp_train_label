use warnings;
use strict;
use Cwd;

my $currentPath = getcwd();# dir for all scripts
my @all_files = `grep -v '^[[:space:]]*\$' $currentPath/QEjobs_status/Dead.txt| grep -v '#'`;#|awk '{print \$2}'`;#all dead QE cases
map { s/^\s+|\s+$//g; } @all_files;
die "No Dead.txt in $currentPath/QEjobs_status" unless(@all_files);
my $submitJobs = "no";
my %sbatch_para = (
            nodes => 1,#how many nodes for your lmp job
            threads => 1,#modify it to 2, 4 if oom problem appears
            #cpus_per_task => 1,#useless if use "mpiexec -np"
            partition => "All",#which partition you want to use
            runPath => "/opt/thermoPW-7-2_intel/bin/pw.x",          
            );

my $jobNo = 1;

for my $i (@all_files){
    print "Job Number $jobNo: $i\n";
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
##SBATCH --cpus-per-task=$sbatch_para{threads}
#SBATCH --partition=$sbatch_para{partition}
##SBATCH --ntasks-per-node=12
##SBATCH --exclude=node23

#source /opt/intel/oneapi/setvars.sh
rm -rf pwscf*
node=$sbatch_para{nodes}
processors=\$(nproc)
threads=\$(nproc)
#threads=$sbatch_para{threads}
np=\$((\$node*\$processors/\$threads))
export OMP_NUM_THREADS=\$threads
#the following two are for AMD CPU if slurm chooses for you!!
export MKL_DEBUG_CPU_TYPE=5
export MKL_CBWR=AUTO
export LD_LIBRARY_PATH=/opt/mpich-4.0.3/lib:/opt/intel/oneapi/mkl/latest/lib:\$LD_LIBRARY_PATH
export PATH=/opt/mpich-4.0.3/bin:\$PATH

/opt/mpich-4.0.3/bin/mpiexec -np \$np $sbatch_para{runPath} -in $basename.in
rm -rf pwscf*
#perl /opt/qe_perl/QEout_analysis.pl
#perl /opt/qe_perl/QEout2data.pl
END_MESSAGE
    unlink "$dirname/$basename.sh";
    open(FH, "> $dirname/$basename.sh") or die $!;
    print FH $here_doc;
    close(FH);
    if($submitJobs eq "yes"){
        chdir($dirname);
        `sbatch $basename.sh`;
        chdir($currentPath);
    }    
}#  

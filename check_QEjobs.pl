=b
Choose which dir you want to check the status of all qe jobs:
my $source_folder = "$currentPath/QE_trimmed4relax";
#my $source_folder = "$currentPath/QEall_set";

=cut

use strict;
use warnings;
use Cwd;
#print "@running\n";
my $whoami = `whoami`;#get username first
$whoami =~ s/^\s+|\s+$//g;
my $currentPath = getcwd();# dir for all scripts
#my $source_folder = "$currentPath/QE_trimmed4relax";#for vc-md

########## source folder you need to assign
my $source_folder = "$currentPath/shear_label/*/labelled";#for vc-relax
#my $source_folder = "$currentPath/thermo_label/*/labelled";#for vc-relax


my @all_QEin = `find $source_folder -type f -name "*.in"`;#keep element info`;
map { s/^\s+|\s+$//g; } @all_QEin;

`rm -rf QEjobs_status`;
`mkdir -p QEjobs_status`;

open(my $FH, "> QEjobs_status/Done.txt") or die $!;
open(my $FH1, "> QEjobs_status/Queueing.txt") or die $!;
open(my $FH2, "> QEjobs_status/Running.txt") or die $!;
open(my $FH3, "> QEjobs_status/Dead.txt") or die $!;

#my @all_paths;
#
#for my $f (@all_QEin){ 
#    my $dir = `dirname $f`;#get path
#    $dir =~ s/^\s+|\s+$//;
#    my $basename = `basename $f`;
#    $basename =~ s/^\s+|\s+$//;
#    $basename =~ s/\.in//g;
#    my $sh_file = "$dir/$basename.sh";
#    push @all_paths, $sh_file;
#}
#
#
##StdOut=/home/jsp/SnPbTe_alloys/dp_train_label/thermo_label/Sn9Pb23Te32-T300-P0-lmpT10to1210-P0-R600000/labelled/lmp_96000.sout
##scontrol show job jobid
#my @jobname_R = `squeue -u $whoami -o "%A %j %u %N %T %M"|grep RUNNING|grep -v JOBID|awk '{print  \$2}'`;#jobnames
my @jobid_R = `squeue -u $whoami -o "%A %j %u %N %T %M"|grep RUNNING|grep -v JOBID|awk '{print  \$1}'`;#jobid
#map { s/^\s+|\s+$//g; } @jobname_R;
map { s/^\s+|\s+$//g; } @jobid_R;

my @running_path;#jobs are currently running
my %path2id;#jobs are currently running
for my $i (@jobid_R){
    my $temp =  `scontrol show job $i|grep StdOut=|awk -F'=' '{print \$2}'`;
    $temp =~ s/^\s+|\s+$//g;
    push @running_path,$temp; 
    $path2id{$temp} = $i;               
}
#            
#my @jobname_Q = `squeue -u $whoami -o "%A %j %u %N %T %M"|grep -v RUNNING|grep -v JOBID|awk '{print  \$2}'`;#jobnames
my @jobid_Q = `squeue -u $whoami -o "%A %j %u %N %T %M"|grep -v RUNNING|grep -v JOBID|awk '{print  \$1}'`;#jobid
#map { s/^\s+|\s+$//g; } @jobname_Q;
map { s/^\s+|\s+$//g; } @jobid_Q;
my @queuing_path;#jobs are currently running
for my $i (@jobid_Q){
    my $temp =  `scontrol show job $i|grep StdOut=|awk -F'=' '{print \$2}'`;
    $temp =~ s/^\s+|\s+$//g;
    push @queuing_path,$temp;                
}
#my @All_jobid = (@jobid_R,@jobid_Q);
#my @all_out_path; #slurm job output path
#for (@all_jobid){
#    my $temp = `scontrol `
#
#}
#
#my @sout_not_in;
#for (@all_paths){
#    print "$_\n";
#
#}
#die;

#print FH $here_doc;
my $doneNu = 0;
my $runNu = 0;
my $queNu = 0;
my $deadNu = 0;

for my $f (@all_QEin){ 
    #my $base = `basename $f`;
    #$base =~ s/^\s+|\s+$//;
    
    #get info from QE input
    my $calculation = `grep calculation $f|awk '{print \$NF}'`;
    $calculation =~ s/^\s+|\s+$//;
    die "No calculation type in $f\n" unless($calculation);

    my $dir = `dirname $f`;#get path
    $dir =~ s/^\s+|\s+$//;
    my $basename = `basename $f`;
    $basename =~ s/^\s+|\s+$//;
    $basename =~ s/\.in//g; 

    my $sh_file = "$dir/$basename.sh";
    #print "$sh_file\n";
    #die;
    open(my $sh, "< $sh_file") or die $!;
    my @tempsh = <$sh>;
    close($sh);
    #get output and job name of a QE job
    my $sout;
    my $jobname;
    #get output filename and job name
    for (@tempsh){
        if(m/#SBATCH\s+--output=\s*(.+)\s*!?/){
            chomp $1;
            $sout = $1;
            die "No QE output name was captured!\n" unless($1);          
        }
        elsif(m/#SBATCH\s+--job-name=\s*(.+)\s*!?/){
            chomp $1;
            $jobname = $1;
            die "No slurm job name was captured!\n" unless($1);          
        }
    }

    if (-e "$dir/$sout"){#sout exists
        my @mark = `grep '!    total energy' $dir/$sout`;
        map { s/^\s+|\s+$//g; } @mark;
        #scf cases
        if($calculation=~m/scf/ and @mark ==1){
            $doneNu++;
            print $FH "$f\n";
            #print "$f\n";
        }
        elsif($calculation=~m/scf/ and @mark != 1){#could be running or dead
            #squeue -o "%A %j %u %N %T %M"
            #398520 jobLi7Al6_mp-1212183-T300-P0 shaohan  PENDING 0:00
            #398523 jobS_mp-77-T50-P0 shaohan node[10,18] RUNNING 1-04:52:12
           # my @jobid_R = `squeue -u $whoami -o "%A %j %u %N %T %M"|grep RUNNING|grep -v JOBID|awk '{print  \$1}'`;#running jobid
           # map { s/^\s+|\s+$//g; } @jobid_R;
            #JobId=412361 JobName=lmp_102000
           

            my $full_path = "$dir/$sout";            
            if($full_path ~~ @running_path){
                my $elapsed = `squeue -u $whoami -o "%A %j %u %N %T %M"|grep $path2id{$full_path}|grep RUNNING`;
                $elapsed =~ s/^\s+|\s+$//g;
               # if($elapsed){
                    $runNu++;
                    print $FH2 "$elapsed for scf:\n$f\n";                    
                    #print "$elapsed for scf:\n$f\n";               
            }
            else{
                $deadNu++;
                print $FH3 "$f\n";
            }
        }

        
    }
    else{#no sout exists, in queue or not submitted
        
        my $full_path = "$dir/$sout";            
        
        if($full_path ~~ @queuing_path){#queneing
            $queNu++;
            print $FH1 "$f\n";
        }
        else{
            $deadNu++;
            print $FH3 "0/0: $f not submitted!\n";#for awk
        }
    }

}
#case fraction
my $total = @all_QEin;
my $doneNuf = $doneNu/$total;
my $runNuf = $runNu/$total;
my $queNuf = $queNu/$total;
my $deadNuf = $deadNu/$total;

print $FH  "#$doneNu/$total: $doneNuf\n";
print $FH2 "#$runNu/$total: $runNuf \n";
print $FH1 "#$queNu/$total: $queNuf \n";
print $FH3 "#$deadNu/$total: $deadNuf\n";

close($FH); 
close($FH1); 
close($FH2); 
close($FH3);
my @deadcases = `cat QEjobs_status/Dead.txt|grep -v '#'`;
map { s/^\s+|\s+$//g; } @deadcases;
if(@deadcases){
    print "\n############\n";
    print "!!!!!The dead cases are:\n";
    system("cat QEjobs_status/Dead.txt");
    print "############\n\n";

}
else{
    print "\n!!!!!No jobs are dead so far!\n\n";
}

my @runningcases = `cat QEjobs_status/Running.txt|grep -v '#'`;
map { s/^\s+|\s+$//g; } @runningcases;
if(@runningcases){
    print "\n++++++++++++\n";
    print "The running cases are:\n";
    system("cat QEjobs_status/Running.txt");
    print "++++++++++++\n";
}
else{
    print "\n!!!!!No jobs are running currently!\n\n";
}

print "\n";
print "***The last line (completed jobs/total jobs) of QEjobs_status/Done.txt!\n";
system("cat QEjobs_status/Done.txt|tail -n 1 ");
print "+++++Check End+++++++\n";

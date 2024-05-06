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
my $source_folder = "$currentPath/thermo_label/*/labelled";#for vc-relax


my @all_QEin = `find $source_folder -type f -name "*.in"`;#keep element info`;
map { s/^\s+|\s+$//g; } @all_QEin;

`rm -rf QEjobs_status`;
`mkdir -p QEjobs_status`;

open(my $FH, "> QEjobs_status/Done.txt") or die $!;
open(my $FH1, "> QEjobs_status/Queueing.txt") or die $!;
open(my $FH2, "> QEjobs_status/Running.txt") or die $!;
open(my $FH3, "> QEjobs_status/Dead.txt") or die $!;

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

    my $nstep = `grep nstep $f|awk '{print \$NF}'`;
    $nstep =~ s/^\s+|\s+$//;
    die "No nstep number in $f\n" unless($nstep);

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
        }
        elsif($calculation=~m/scf/ and @mark != 1){
            #squeue -o "%A %j %u %N %T %M"
            #398520 jobLi7Al6_mp-1212183-T300-P0 shaohan  PENDING 0:00
            #398523 jobS_mp-77-T50-P0 shaohan node[10,18] RUNNING 1-04:52:12
            my @submitted = `squeue -u $whoami -o "%A %j %u %N %T %M"|awk '{print  \$2}'`;#jobnames
            my @submitted1 = `squeue -u $whoami -o "%A %j %u %N %T %M"|awk '{print  \$1}'`;#jobid
            map { s/^\s+|\s+$//g; } @submitted;
            map { s/^\s+|\s+$//g; } @submitted1;
            my %jobname2id;
            @jobname2id{@submitted} = @submitted1;            
            #for my $t (@submitted1){
            #    $t =~ m/(\d+)\s+(.+)/;
            #    chomp ($1,$2);
            #    $jobname2id{$2} = $1;
            #}
            if($jobname ~~ @submitted){
                my $elapsed = `squeue|grep $jobname2id{$jobname}`;
                $elapsed =~ s/^\s+|\s+$//g;
               # if($elapsed){
                    $runNu++;
                    print $FH2 "$elapsed for scf:\n$f\n";
               # }
               # else{
               #     $deadNu++;
               #     print $FH3 "$f\n";
               # }
            }
            else{
                $deadNu++;
                print $FH3 "$f\n";
            }
        }

        #md cases
        if($calculation=~m/(md|relax)/ and @mark == $nstep){
            $doneNu++;
            print $FH "$f\n";
        }
        elsif($calculation=~m/(md|relax)/ and @mark < $nstep){
            #squeue -o "%A %j %u %N %T %M"
            #398520 jobLi7Al6_mp-1212183-T300-P0 shaohan  PENDING 0:00
            #398523 jobS_mp-77-T50-P0 shaohan node[10,18] RUNNING 1-04:52:12
            my @submitted = `squeue -u $whoami -o "%A %j %u %N %T %M"|awk '{print  \$2}'`;#jobnames
            my @submitted1 = `squeue -u $whoami -o "%A %j %u %N %T %M"|awk '{print  \$1}'`;#jobid
            map { s/^\s+|\s+$//g; } @submitted;
            map { s/^\s+|\s+$//g; } @submitted1;
            my %jobname2id;
            @jobname2id{@submitted} = @submitted1;

            if($jobname ~~ @submitted){#running
                my $elapsed = `squeue|grep $jobname2id{$jobname}`;
                #if($elapsed){
                    $elapsed =~ s/^\s+|\s+$//g;                
                    $runNu++;
                    my $temp = @mark."/".$nstep;
                    print $FH2 "**$elapsed\n $temp: in $f\n\n";
                #}
                #else{
                #    $deadNu++;
                #    my $temp = @mark."/".$nstep;
                #    print $FH3 "$temp: $f !\n";#for awk    
                #}
            }
            else{
                $deadNu++;
                my $temp = @mark."/".$nstep;
                print $FH3 "$temp: $f !\n";#for awk
            }
        }
    }
    else{#no sout exists, in queue
         my @submitted = `squeue -u $whoami -o "%A %j %u %N %T %M"|awk '{print  \$2}'`;#jobnames
            map { s/^\s+|\s+$//g; } @submitted;
        if($jobname ~~ @submitted){#queneing
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



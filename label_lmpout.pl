=b

=cut
use warnings;
use strict;
use JSON::PP;
use Data::Dumper;
use List::Util qw(min max any shuffle);
use Cwd;
use POSIX;
use Parallel::ForkManager;
#use lib '../scripts';#assign pm dir
#use elements;#all setting package

my $currentPath = getcwd();# dir for all scripts
chdir("..");
my $mainPath = getcwd();# main path of Perl4dpgen dir
chdir("$currentPath");

###settings here!
my $max4relabel = 100;# how many cfgs you want to use in a labelled folder
my $lowerbound = 0.03;#below which not lebelled
my $upperbound = 0.3;#above which not lebelled
my $sour_folder = "$currentPath/shear_label";
#my $sour_folder = "$currentPath/thermo_label";

my $forkNo = 1;#although we don't have so many cores, only for submitting jobs into slurm
my $pm = Parallel::ForkManager->new("$forkNo");

my @all_subfolders = `find $sour_folder -maxdepth 1 -mindepth 1 -type d -name "*" -exec readlink -f {} \\;|sort`;
map { s/^\s+|\s+$//g; } @all_subfolders;
my $sub_No = @all_subfolders;
my $ng_count = 0;
`rm -f $sour_folder/nolabel.dat`;
`touch $sour_folder/nolabel.dat`;
for my $n (@all_subfolders){#loop over folders with no labelled sub folder
    `rm -rf $n/labelled`;
    
    unless (-e "$n/md.out"){
        #print "no md.out in $n. Skipped!\n";
       `echo "no md.out in $n" >> $sour_folder/nolabel.dat`;
       $ng_count++;
        next; 
    }
    #  step   max_devi_v   min_devi_v   avg_devi_v   max_devi_f  min_devi_f         avg_devi_f
    my @maxf = `grep -v '^[[:space:]]*\$' $n/md.out | grep -v step|awk '{print \$5}'`;
    map { s/^\s+|\s+$//g; } @maxf;
    my @step_maxf1 = `grep -v '^[[:space:]]*\$' $n/md.out | grep -v step|awk '{print \$1 " " \$5}'`;
    map { s/^\s+|\s+$//g; } @step_maxf1;
    my @step_maxf = @step_maxf1;
    #my @step_maxf = shuffle @step_maxf1;
    unless (@maxf){#
        print "No max force deviation in $n\n";
       `echo "No max force deviation in $n" >> $sour_folder/nolabel.dat`;
       $ng_count++;
       next; 
    }
    #max forces within criterions
    my @force_dev = grep { $_ >= $lowerbound and $_ <= $upperbound}  @maxf; 
    unless (@force_dev){#
        print "No max force deviation between $lowerbound and $upperbound in $n\n";
        `echo "No max force deviation between $lowerbound and $upperbound in $n" >> $sour_folder/nolabel.dat`;
        $ng_count++;
        next; 
    }

    #my @step = `grep -v '^[[:space:]]*\$' $n/md.out | grep -v step|awk '{print \$1}'`;
    #map { s/^\s+|\s+$//g; } @step;        
    #my @maxsort = sort {$a <=> $b} @maxf;
   # my $max4relabel
    
    `mkdir $n/labelled`;
    ####begin relabelling
    my $counter = 0;
    for my $m (0..$#step_maxf){#original sequence
    #for my $m (0..$#maxf){#original sequence
            #print "## $m: $step_maxf[$m]\n";
        my @temp = split(/\s+/,$step_maxf[$m]);
        #if(any { $_ == $maxf[$m] } @force_dev){    
        if(any { $_ == $temp[1] } @force_dev and $counter < $max4relabel){    
           $counter++;
           my $lmpfile = "lmp_$temp[0].cfg";
           `cp $n/lmp_output/$lmpfile $n/labelled/$lmpfile`
        }
    }
}

my @all_labelled = `find $sour_folder/*/labelled -maxdepth 1 -mindepth 1 -type f -name "*.cfg"`;
map { s/^\s+|\s+$//g; } @all_labelled;
my $labelNo = @all_labelled;

print "Total labelled cfg Number: $labelNo\n";
print "If the above number is too large for SCF, you can reduce \$upperbound to 0.20 or 0.18...\n";
#for (@all_labelled){
#    print "$_\n";
#}
#die;
`echo " " >> $sour_folder/nolabel.dat`;
`echo "### label process summary:" >> $sour_folder/nolabel.dat`;
`echo "folders without labelled subfolder/total folders => $ng_count/$sub_No" >> $sour_folder/nolabel.dat`;
system("cat $sour_folder/nolabel.dat");

print "\nFor details, check $sour_folder/nolabel.dat!!!\n";

#else{#folders without doing label process before (first time to do label)
#    print "###no nolabel.dat\n\n";
#    #get all subfolder names under source folder
#    my @temp = `find $sour_folder -maxdepth 1 -mindepth 1 -type d -name "*" -exec readlink -f {} \\;|sort`;
#    map { s/^\s+|\s+$//g; } @temp;
#    my @nolabel;
#    my $count = 0;
#    for my $d (@temp){
#        unless(-e "$d/labelled"){
#            $count++;
#            `touch $sour_folder/nolabel.dat` if($count == 1);
#            push @nolabel,$d;
#            `echo $d >> $sour_folder/nolabel.dat`;
#        }
#    }
#    my $nolabelNo = @nolabel;
#    unless($nolabelNo){
#        print "all folders with labelled subfolders. No need to relabel.\n";
#    }
#    else{
#        print "##Folders without labelled subfolders: $nolabelNo\n";
#        print "\n\n!!!!!! PLEASE CHECK nolabel.dat and then conduct the same script again after setting a reasonable larger \$upperbound (<0.3).\n";
#    }
#}
#
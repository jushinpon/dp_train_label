=b

=cut
use warnings;
use strict;
use JSON::PP;
use Data::Dumper;
use List::Util qw(min max any);
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
my $sour_folder = "$currentPath/thermo_label";

my @all_subfolders = `find $sour_folder -maxdepth 1 -mindepth 1 -type d -name "*" -exec readlink -f {} \\;|sort`;
map { s/^\s+|\s+$//g; } @all_subfolders;
my $sub_No = @all_subfolders;
`rm -f $sour_folder/no_mdout.dat`;
`touch $sour_folder/no_mdout.dat`;
for my $n (@all_subfolders){#loop over folders with no labelled sub folder
    unless (-e "$n/md.out"){
        
        print "no md.out in $n. Skipped!\n";
       `echo "$n" >> $sour_folder/no_mdout.dat`;
    }
}
print "printing no_mdout.dat!\n\n";
system("cat $sour_folder/no_mdout.dat");

print "\nFor details, check $sour_folder/no_mdout.dat!!!\n";

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
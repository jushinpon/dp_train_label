=b

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

######### modify the sou_dir to the correct one (thermo_label, sfe_label, surf_label)
my $no_mdout = "$currentPath/thermo_label/no_mdout.dat";#source file

my @all_subfolders = `cat $no_mdout|grep -v '^[[:space:]]*\$'`;
map { s/^\s+|\s+$//g; } @all_subfolders;

for my $i (@all_subfolders){
    my $basename = `basename $i`;
    $basename =~ s/^\s+|\s+$//g;

    my $dirname = `dirname $i`;
    $dirname =~ s/^\s+|\s+$//g;
   # print "$dirname/$basename\n";
    #system("ls $dirname/$basename/$basename.sh");
    chdir("$dirname/$basename");
    `sbatch $basename.sh`;        
}#  


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
my $sou_dir = "$currentPath/thermo_label";#source dir

my $forkNo = 1;#although we don't have so many cores, only for submitting jobs into slurm
my $pm = Parallel::ForkManager->new("$forkNo");

my @all_sh = `find $sou_dir -type f -name "*.sh" -exec readlink -f {} \\;|sort`;
map { s/^\s+|\s+$//g; } @all_sh;

for my $i (@all_sh){
    my $basename = `basename $i`;
    $basename =~ s/^\s+|\s+$//g;

    my $dirname = `dirname $i`;
    $dirname =~ s/^\s+|\s+$//g;
    chdir($dirname);
    `sbatch $basename`;        
}#  


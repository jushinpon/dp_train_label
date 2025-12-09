=b
choose only one of the following two:

#y $filefold = "QE_trimmed4relax";
my $filefold = "QEall_set";
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
#my $sou_dir = "$currentPath/scale_label/*/labelled";#source dir with labelled subfolders
#my $sou_dir = "$currentPath/shear_label/*/labelled";#source dir with labelled subfolders
my $sou_dir = "$currentPath/thermo_label/*/labelled";#source dir with labelled subfolders

my $forkNo = 1;#although we don't have so many cores, only for submitting jobs into slurm
my $pm = Parallel::ForkManager->new("$forkNo");

my @all_files = `find $sou_dir -maxdepth 1 -mindepth 1 -type f -name "*.sh" `;
#my @all_files = `find $currentPath/$filefold -maxdepth 2 -mindepth 2 -type f -name "*.sh" -exec readlink -f {} \\;|sort`;
map { s/^\s+|\s+$//g; } @all_files;

for my $i (@all_files){
    my $dirname = `dirname $i`;
    $dirname =~ s/^\s+|\s+$//g;
    chdir($dirname);
    my $prefix = `basename $i`;
    $prefix =~ s/^\s+|\s+$//g;
    $prefix =~ s/\.sh//g;
    unlink "$prefix.sout";
    `sbatch $prefix.sh`;
}#  


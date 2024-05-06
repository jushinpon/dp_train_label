=b
make lmp input files for all strucutres in labelled folders.
You need to use this script in the dir with all dpgen collections (in all_cfgs folder)
perl ../tool_scripts/cfg2lmpinput.pl 
=cut
use warnings;
use strict;
use Cwd;
use POSIX;

my $currentPath = getcwd();# dir for all scripts
my @allQEin = `grep -v '^[[:space:]]*\$' $currentPath/QEjobs_status/Dead.txt| grep -v '#'|awk '{print \$2}'`;#all dead QE cases
map { s/^\s+|\s+$//g; } @allQEin;

for my $i (@allQEin){
    #print "\$i: $i\n";

    my $dirname = `dirname $i`;
    $dirname =~ s/^\s+|\s+$//g;
    chdir($dirname);
    my $prefix = `basename $i`;
    $prefix =~ s/^\s+|\s+$//g;
    $prefix =~ s/\.in//g;
    unlink "$prefix.sout";
    `sbatch $prefix.sh`;
}#  


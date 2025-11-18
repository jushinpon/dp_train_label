#!/usr/bin/perl
=b
=cut
use strict;
use Cwd;
use Data::Dumper;
use JSON::PP;
use Data::Dumper;
use List::Util qw(min max);
use Cwd;
use POSIX;
use Parallel::ForkManager;

# You need to decide what kind of label process you want to do 
#1. thermodynamics 
#2. surface
#3. stacking fault energy (SFE)


#my @label_types = ("data4thermo");#,"surface","SFE");
# remove old folders and 
`rm -rf data4shear`;
`mkdir -p data4shear`;

###parameters to set first
my $currentPath = getcwd();# dir for all scripts
chdir("..");
my $mainPath = getcwd();# main path of Perl4dpgen dir
chdir("$currentPath");

#######pick what you want for the label process
my @datafiles = `find -L ../dp_train_new/initial -type f -name "*.data"|grep Sn12Pb20Te32-T300-P0.data`;#|grep -v 111|grep -v 110|grep -v 100`;#find all data files
map { s/^\s+|\s+$//g; } @datafiles;
die "No data files to collect!\n" unless(@datafiles);

for (@datafiles){
    print "$_\n";
    `cp $_ data4shear/`;
}

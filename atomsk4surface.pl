#!/usr/bin/perl

#"[010] [001] [100]";# 100 crystal axis vectors,
#"[001] [1-10] [110]";# 110 crystal axis vectors,x 2
#"[1-21] [10-1] [111]";# 111 crystal axis vectors, x6
use strict;
use warnings;
use Cwd;

my $currentPath = getcwd();# dir for all scripts
chdir("..");
my $mainPath = getcwd();# main path of Perl4dpgen dir
chdir("$currentPath");

my %scale = (
    111 => 5,#scale 1/5 z length of cell
    110 => 2,# x 2
    100 => 1
);

my $rcut = 6; #Your DLP rcut
my $source = "$currentPath/data4surface";
my $des = "$currentPath/surface_data";

my @surfaceZ = (
    '[010] [001] [100]',
    '[001] [1-10] [110]',
    '[1-21] [10-1] [111]'
);
my @sur_name = map {my @temp = split (/\s+/,$_); my $last = $temp[-1]; $last =~ s/[\[\]]//g; $last} @surfaceZ;

`rm -rf $des`;
`mkdir -p $des`;

my @source_data = `find -L $source -mindepth 1 -maxdepth 1 -type f -name "*.data" `;
map { s/^\s+|\s+$//g; } @source_data;

for my $file (@source_data) {
    my $basename = `basename $file`;
    my $dirname = `dirname $file`;
    $basename =~ s/^\s+|\s+$//g;
    $dirname =~ s/^\s+|\s+$//g;
    $basename =~ s/\.data//g;
    `rm -f $dirname/$basename.lmp`;#for atomsk only
    `cp $dirname/$basename.data $dirname/$basename.lmp`;#for atomsk only

    for my $i (0..$#sur_name){
        `rm -f $dirname/$basename-$sur_name[$i].lmp`;
        my $command = "atomsk $dirname/$basename.lmp -orient [100] [010] [001] $surfaceZ[$i] -orthogonal-cell -reduce-cell $des/$basename-$sur_name[$i].lmp";
        system("$command");
        my @x = `grep "xlo xhi" $des/$basename-$sur_name[$i].lmp|awk '{print \$1 "\\n" \$2}'`;
        map { s/^\s+|\s+$//g; } @x;
        
        my @y = `grep "ylo yhi" $des/$basename-$sur_name[$i].lmp|awk '{print \$1 "\\n" \$2}'`;
        map { s/^\s+|\s+$//g; } @y;
        
        my @z = `grep "zlo zhi" $des/$basename-$sur_name[$i].lmp|awk '{print \$1 "\\n" \$2}'`;
        map { s/^\s+|\s+$//g; } @z;
        my $zleng = $z[1] - $z[0];
        
        my $ref_hi = $z[0] + $zleng * (1./$scale{$sur_name[$i]});

        system("atomsk $des/$basename-$sur_name[$i].lmp -cut above $ref_hi z $des/$basename-$sur_name[$i]_cut.lmp");
        `rm -f $des/$basename-$sur_name[$i].lmp`;

        my $final_zhi = $ref_hi + $rcut;
        my $box = "$final_zhi z ";
        #my $box = "$x[0] $x[1] $y[0] $y[1] $z[0] $final_zhi ";
        `sed -i '/zlo zhi/c\  $z[0]  $final_zhi zlo zhi' $des/$basename-$sur_name[$i]_cut.lmp`;
        `mv $des/$basename-$sur_name[$i]_cut.lmp $des/$basename-$sur_name[$i].data`;
        `rm -f $des/$basename-$sur_name[$i]_cut.lmp`;

    }
    `rm -f $des/$basename.lmp`;#remove old
}


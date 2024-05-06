=b
make qe input files for all strucutres in labelled folders.
You need to use this script in the dir with all dpgen collections (in all_cfgs folder)
perl ../tool_scripts/cfg2QEinput.pl 
=cut
use warnings;
use strict;
use JSON::PP;
use Data::Dumper;
use List::Util qw(min max);
use Cwd;
use POSIX;
use Parallel::ForkManager;
use lib '.';#assign pm dir
use elements;#all setting package
use List::Util qw/shuffle/;

my $currentPath = getcwd();# dir for all scripts
chdir("..");
my $mainPath = getcwd();# main path of Perl4dpgen dir
chdir("$currentPath");

########## source folder you need to assign
my $sou_dir = "$currentPath/thermo_label";#source dir with labelled subfolders
my $QEin_template_dir = "$mainPath/dp_train_new/initial";#the folders to find the QE input file for using the same setting.
my $data_dir = "$currentPath/data4thermo";#the prefix of data file is used to find the corresponding folder in $QEin_template_dir.

my $forkNo = 1;#although we don't have so many cores, only for submitting jobs into slurm
my $pm = Parallel::ForkManager->new("$forkNo");

#you may use relabel script to increase more labelled folders if needed
my @labelled_dir = `find $sou_dir -type d -name "labelled" -exec readlink -f {} \\;|sort`;
map { s/^\s+|\s+$//g; } @labelled_dir;
my @pathOfAllcfgs;
for my $dir (@labelled_dir){
    #print "$dir\n";
    my @cfgs = <$dir/*.cfg>;
    map { s/^\s+|\s+$//g; } @cfgs;
    for my $c (sort @cfgs){
        push @pathOfAllcfgs,$c;
    }
}#  

my $labelledNo = @pathOfAllcfgs;
print "***all labelled cfg Number: $labelledNo.\n";

#make qe input files
my $counter = 0;
for my $cfg (@pathOfAllcfgs){
    $counter++;
    my $basename = `basename $cfg`;
    $basename =~ s/^\s+|\s+$//g; 

    my $dirname = `dirname $cfg`;
    $dirname =~ s/^\s+|\s+$//g; 
    $basename =~ s/\.cfg//g; 
    `rm -f $dirname/$basename.in`;#remove the old QE input
    #print "\$basename,\$dirname,\$folder: $basename,$dirname\n";
    #print "No $counter: $cfg\n";
    #find reference QE input in /initial folder
    $cfg =~ m#.+/(.+T\d+-P\d+)-lmp.+/labelled/lmp_\d+.cfg#;    
    #print "\$1: $1\n";

    my @qein = `grep -v '^[[:space:]]*\$' $QEin_template_dir/$1/$1.in`;
    die "no qe input for $QEin_template_dir/$1/$1.in\n" unless(@qein);
    map { s/^\s+|\s+$//g; } @qein;
    #find line number of key words
    my $cal_line;
    my $cell_line;
    my $pos_line;
    my $nat;
    for my $q (0..$#qein){
        my $temp = $qein[$q];
        if($temp =~ /^calculation/){$cal_line = $q;}
        elsif($temp =~ /^ATOMIC_POSITIONS/){$pos_line = $q;}#print "$q\n";
        elsif($temp =~ /^CELL_PARAMETERS/){$cell_line = $q;}#print "$q\n";
        elsif($temp =~ /^nat\s*=\s*(\d*)/){chomp $1;$nat = $1;}#print "$q\n";
    }

    #modify calculation type
    $qein[$cal_line] = "calculation = \'scf\'";
    #get cell information from cfg
    my $xlo_bound = `sed -n '/ITEM: BOX BOUNDS.*/{n;p}' $cfg | awk '{print \$1}'`;
    my $xhi_bound = `sed -n '/ITEM: BOX BOUNDS.*/{n;p}' $cfg | awk '{print \$2}'`;
    my $xy = `sed -n '/ITEM: BOX BOUNDS.*/{n;p}' $cfg | awk '{print \$3}'`;
    chomp ($xlo_bound,$xhi_bound,$xy);
    unless ($xy){$xy = 0.0;}#if no $xy
    my $ylo_bound = `sed -n '/ITEM: BOX BOUNDS.*/{n;n;p}' $cfg | awk '{print \$1}'`;
    my $yhi_bound = `sed -n '/ITEM: BOX BOUNDS.*/{n;n;p}' $cfg | awk '{print \$2}'`;
    my $xz = `sed -n '/ITEM: BOX BOUNDS.*/{n;n;p}' $cfg | awk '{print \$3}'`;
    chomp ($ylo_bound,$yhi_bound,$xz);
    unless ($xz){$xz = 0.0;}#if no $xz
    my $zlo_bound = `sed -n '/ITEM: BOX BOUNDS.*/{n;n;n;p}' $cfg | awk '{print \$1}'`;
    my $zhi_bound = `sed -n '/ITEM: BOX BOUNDS.*/{n;n;n;p}' $cfg | awk '{print \$2}'`;
    my $yz = `sed -n '/ITEM: BOX BOUNDS.*/{n;n;n;p}' $cfg | awk '{print \$3}'`;
    chomp ($zlo_bound,$zhi_bound,$yz);
    unless ($yz){$yz = 0.0;}#if no $yz
    my $xlo = $xlo_bound - min(0.0,$xy,$xz,($xy+$xz));
    my $xhi = $xhi_bound - max(0.0,$xy,$xz,($xy+$xz));
    my $ylo = $ylo_bound - min(0.0,$yz);
    my $yhi = $yhi_bound - max(0.0,$yz);
    my $lx = sprintf("%.6f",$xhi - $xlo) ;
    my $ly = sprintf("%.6f",$yhi - $ylo) ;
    my $lz = sprintf("%.6f",$zhi_bound - $zlo_bound);
    my $zero = sprintf("%.6f",0.0);
    $xy = sprintf("%.6f",$xy); 
    $xz = sprintf("%.6f",$xz); 
    $yz = sprintf("%.6f",$yz); 
    #modify qe cell below
    $qein[$cell_line + 1] =  "$lx $zero $zero";
    $qein[$cell_line + 2] =  "$xy $ly $zero";
    $qein[$cell_line + 3] =  "$xz $yz $lz";

    #coordinates of cfg
    my @coor_cfg = `grep -v '^[[:space:]]*\$' $cfg|grep -A $nat "ITEM: ATOMS"|grep -v "ITEM: ATOMS"|grep -v -- '--'|awk '{print \$(NF-2)\" \" \$(NF-1)\" \" \$NF}'`;
    map { s/^\s+|\s+$//g; } @coor_cfg; 
    for my $l ($pos_line+1 ..$pos_line+$nat){
        my @temp = split(/\s+/,$qein[$l]);
        chomp $temp[0];#element
        #print "***$l: $qein[$l]\n";
        $qein[$l] = "$temp[0] ". "$coor_cfg[$l - $pos_line - 1]";
        #print "$l: $qein[$l]\n";
    }
    my $qeinput = join("\n",@qein);
    chomp $qeinput;

    open(FH, "> $dirname/$basename.in") or die $!;
    print FH $qeinput;
    close(FH);    
}

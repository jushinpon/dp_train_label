#!/usr/bin/perl
=b
usage: perl QEoutput2data.pl [sout dir]

if no sout dir, the current dir will be used to find sout files.
=cut

#density (g/cm3), arrangement, mass, lat a , lat c
    #my @temp = 
#    @{$used_element{$_}} = &elements::eleObj("$_");
use strict;
use warnings;
use Data::Dumper;
use POSIX;
use lib '.';
use elements;#all setting package
use Cwd;


my $currentPath = getcwd();
my @QEin_files = `find -L ../initial -type f -name "*.in"`;#sout in current folder
map { s/^\s+|\s+$//g; } @QEin_files;
unless(@QEin_files){die "no QE input files to convert into lmp data files.\n";}
print "total folders after $0: ".@QEin_files."\n";

for my $file (@QEin_files){
    my $QEin_path = `dirname $file`;
    $QEin_path =~ s/^\s+|\s+$//g;
    my $QEin_name = `basename $file`;
    $QEin_name =~ s/^\s+|\s+$//g;
    $QEin_name =~ s/\..*//g;
    my $natom = `grep "nat" $file|awk '{print \$3}'`;
    die "No atom number was found in $file" unless ($natom); 
    $natom =~ s/^\s+|\s+$//g;
    my $ntype = `grep "ntyp" $file|awk '{print \$3}'`;
    die "No atom type was found in $file" unless ($ntype); 
    $ntype =~ s/^\s+|\s+$//g;
    #get cell information of all frames
    #CELL_PARAMETERS (angstrom) or {angstrom}
    #5.631780735   0.001261244   0.001887268
    my @CELL_PARAMETERS = `grep -A3 "CELL_PARAMETERS" $file|grep -v "CELL_PARAMETERS"|grep -v -- '--'`;
    map { s/^\s+|\s+$//g; } @CELL_PARAMETERS; 
    die "No CELL_PARAMETERS were found in $file" unless (@CELL_PARAMETERS); 
    #my @box;#array of array, equal to frame numbers
    #get atom coords information of all frames
    #ATOMIC_POSITIONS (angstrom)
    #Co            2.7414458575        2.7928470261        2.8314219861
    my @coords = `grep -A $natom "ATOMIC_POSITIONS" $file|grep -v "ATOMIC_POSITIONS"|grep -v -- '--'`;
    map { s/^\s+|\s+$//g; } @coords;
    die "coord info is not available in $file\n" unless(@coords);
    #element types of atoms for all frames 
    my @Alltypes = `grep -A $natom "ATOMIC_POSITIONS" $file|grep -v "ATOMIC_POSITIONS"|awk '{print \$1}'|grep -v -- '--'`;
    map { s/^\s+|\s+$//g; } @Alltypes;
    my @element4atoms = @Alltypes;#only need the first set information
    my @used_ele = `grep -A $ntype "ATOMIC_SPECIES" $file|grep -v "ATOMIC_SPECIES"|awk '{print \$1}'|grep -v -- '--'`; 
    map { s/^\s+|\s+$//g; } @used_ele;    
    my %ele2id = map { $used_ele[$_] => $_ + 1  } 0 .. $#used_ele;#make a hash for element -> type id for lmp
    #1 65.409  # Co
    #2 15.9994  # Ni
    my $mass4data;#mass for each element
    for my $t (0..$#used_ele){
        my $ele = $used_ele[$t];
        #density (g/cm3), arrangement, mass        
        my $mass = &elements::eleObj("$ele")->[2];
        $mass4data .= $t+1 . " $mass  \# $ele\n";
    }
    
    chomp $mass4data;#move the new line for the last line
    # making data files below:

    my @box;
    for my $d (0..2){#loop over cell vectors
        my @temp = split (/\s+/,$CELL_PARAMETERS[$d]);
        map { s/^\s+|\s+$//g; } @temp;
        for my $i (@temp){ #loop over each component of each vector
            push @{$box[$d]},$i;
        }
    }  
    
    my $a = ( ${$box[0]}[0]**2 + ${$box[0]}[1]**2 + ${$box[0]}[2]**2 )**0.5;
    my $b = ( ${$box[1]}[0]**2 + ${$box[1]}[1]**2 + ${$box[1]}[2]**2 )**0.5;
    my $c = ( ${$box[2]}[0]**2 + ${$box[2]}[1]**2 + ${$box[2]}[2]**2 )**0.5;
    my $cosalpha = (${$box[1]}[0]*${$box[2]}[0] + ${$box[1]}[1]*${$box[2]}[1] + ${$box[1]}[2]*${$box[2]}[2])/($b*$c);
    my $cosbeta  = (${$box[2]}[0]*${$box[0]}[0] + ${$box[2]}[1]*${$box[0]}[1] + ${$box[2]}[2]*${$box[0]}[2])/($c*$a);
    my $cosgamma = (${$box[0]}[0]*${$box[1]}[0] + ${$box[0]}[1]*${$box[1]}[1] + ${$box[0]}[2]*${$box[1]}[2])/($a*$b);
    my $lx = sprintf("%.6f",$a);
    my $xy = sprintf("%.6f",$b*$cosgamma);
    my $xz = sprintf("%.6f",$c*$cosbeta);
    my $ly = sprintf("%.6f",sqrt($b**2 - $xy**2));
    my $yz = sprintf("%.6f",($b*$c*$cosalpha-$xy*$xz)/$ly);
    my $lz = sprintf("%.6f",sqrt($c**2 - $xz**2 - $yz**2));
    my @cell4data = (
        "0.000000 $lx xlo xhi",
        "0.000000 $ly ylo yhi",
        "0.000000 $lz zlo zhi",
        "$xy $xz $yz xy xz yz"
    );
    map { s/^\s+|\s+$//g; } @cell4data;     
    my $cell4data = join("\n",@cell4data);
    my $coords4data;
    #print "$cell4data\n";
    for my $d (0..$natom -1){#loop over coords
         my @tempQEcoord = split (/\s+/,$coords[$d]);
         map { s/^\s+|\s+$//g; } @tempQEcoord;
         #print "@tempQEcoord[1..3]\n";
        # for my $i (@temp){ #loop over each component of each vector
             my $temp_coord = eval($d + 1)." $ele2id{$tempQEcoord[0]} @tempQEcoord[1..3]\n";
             #my $temp_coord = ($d + 1) . " "."$ele2id{$tempQEcoord[0]}" . " ". "@tempQEcoord[1..3]\n";
             $coords4data .= $temp_coord;
         # }
        }#over QE atom coords
        chomp $coords4data;
### set hash for heredoc
        my %hash_para = (
        output_file => " $QEin_path/$QEin_name.data",
        natom => "$natom",
        ntype => "$ntype",
        cell => "$cell4data",
        masses => "$mass4data",            
        coords => "$coords4data"
        );
        &make_data_file(\%hash_para);
}#all QE input files

sub make_data_file{

my ($para_hr) = @_;
my $here_doc =<<"END_MESSAGE";
# LAMMPS data file written by OVITO Basic 3.7.8

$para_hr->{natom} atoms
$para_hr->{ntype} atom types

$para_hr->{cell}

Masses

$para_hr->{masses}

Atoms  # atomic

$para_hr->{coords}
END_MESSAGE
my $temp = $para_hr->{output_file};
chomp $temp;
open(FH, "> $temp") or die $!;
print FH $here_doc;
close(FH);
#`cat << $QEinput > $temp`;
}

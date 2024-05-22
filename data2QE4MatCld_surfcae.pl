#!/usr/bin/perl
=b

If you want to use your own QE setting, you need to modify the heredoc.

https://docs.lammps.org/Howto_triclinic.html
The parallelepiped has its “origin” at (xlo,ylo,zlo) and 
is defined by 3 edge vectors starting from the origin given by
a = (xhi-xlo,0,0); b = (xy,yhi-ylo,0); c = (xz,yz,zhi-zlo). 

cfg format:
ITEM: BOX BOUNDS xy xz yz
xlo_bound xhi_bound xy
ylo_bound yhi_bound xz
zlo_bound zhi_bound yz

data format:
0.000000000000      11.327003370000  xlo xhi
0.000000000000      10.970021530071  ylo yhi
0.000000000000      10.529870936067  zlo zhi
0.507912903230       0.714988108603      -0.046041891158  xy xz yz
=cut
use strict;
use warnings;
use Cwd;
use Data::Dumper;
use JSON::PP;
use Data::Dumper;
use List::Util qw(min max);
use Cwd;
use POSIX;
use Parallel::ForkManager;
use lib './';#assign pm dir
use elements;#all setting package

print "\n\n***IMPORTANT!!! If all your data files do not include all elements you want to use\n";
print ", you need to modify heredoc for setting the largest ecutrho and ecutwfc of your element. \n";
print "\nMaybe you also need to modify cell_dofree setting for your QE cases. \n\n";
###parameters to set first
my $currentPath = getcwd();
#my @myelement =  ("B","N");#corresponding to lmp type ids
my @datafile = `find $currentPath/surface_data -name "*.data"`;#find all data files
chomp @datafile;
die "No data files\n" unless(@datafile);
#collect all elements in data files
my %myelement;#get all ements from all data files
for (@datafile){
    #1  atom types
    my $typeNum = `grep "atom types" $_|awk '{print \$1}'`;
    chomp  $typeNum;
    #1   58.93319500             # Co
    my @ele = `grep -v '^[[:space:]]*\$' $_|grep -A $typeNum Masses|grep -v Masses|grep -v -- '--'|awk '{print \$NF}'`;
    #my @ele = `grep -v '^[[:space:]]*\$' $_|grep -A $typeNum Masses|grep -v Masses|grep -v -- '--'|awk '{print \$NF}'`;
    map { s/^\s+|\s+$//g; } @ele;
    die "No Masses for finding element symbol in $_\n" unless(@ele);
    @myelement{@ele} = 1;
}

my @myelement = keys  %myelement;
map { s/^\s+|\s+$//g; } @myelement;
die "No elements were found\n" unless (@myelement);
my @temperature = ("10");#temperatures for QE_MD, only template for the following sed trim
my @pressure = ("0");#pressure for vc-md, only template for the following sed trim
my $calculation = "vc-md";#set temperature and pressure to be 0 0 for scf
my $stepsize = 20;#20 ~ 0.97 fs
my $nstep = 300;#how many steps for md for vc-relax
my $pseudo_dir = "/opt/QEpot/SSSP_efficiency_pseudos/";
####end of setting parameters
###get pot setting here!
my $json;
{
    local $/ = undef;
    open my $fh, '<', '/opt/QEpot/SSSP_efficiency.json' or die "no QE pot file\n";#or precision
    $json = <$fh>;
    close $fh;
}
my $decoded = decode_json($json);
#print "$decoded->{\"B\"}->{filename}\n";
#die;
######   cutoff   #######
my @rho_cutoff;
my @cutoff;
my %pot;
#my @pot;
my %used_element;
for (@myelement){#unique
    chomp;
    #print "element: $_\n";
    push @rho_cutoff,$decoded->{$_}->{rho_cutoff};
    push @cutoff,$decoded->{$_}->{cutoff};
 #density (g/cm3), arrangement, mass, lat a , lat c
    #my @temp = 
    $used_element{$_} = &elements::eleObj("$_");#return reference only
    #@{$used_element{$_}} = &elements::eleObj("$_");#return reference only
    die "no information of element $_ in elements.pm\n" unless ($used_element{$_});
    die "no pot file of element $_\n" unless ($decoded->{"$_"}->{filename});
    my $temp1 = ${$used_element{$_}}[2];
    my $temp2 = $decoded->{$_}->{filename};
    my @temp = ($_,$temp1,$temp2);
    $pot{$_} = join(" ",@temp);
    chomp $pot{$_};
}
# for keeping the largest ones only
my @rho_cutoff1 = sort {$b<=>$a} @rho_cutoff;
my @cutoff1 = sort {$b<=>$a} @cutoff;
#my $ntyp = @pot; #element type number
## for here doc
my $rho_cutoff = $rho_cutoff1[0];
my $cutoff = $cutoff1[0];

`rm -rf data2QE4MatCld`;
`mkdir data2QE4MatCld`;

for my $id (@datafile){
    #get required data first
    my %current_elements;#get all ements of a data file
    my $typeNum = `grep "atom types" $id|awk '{print \$1}'`;
    chomp  $typeNum;
    my @ele = `grep -v '^[[:space:]]*\$' $id|grep -A $typeNum Masses|grep -v Masses|grep -v -- '--'|awk '{print \$NF}'`;
    #my @ele = `grep -v '^[[:space:]]*\$' $_|grep -A $typeNum Masses|grep -v Masses|grep -v -- '--'|awk '{print \$NF}'`;
    map { s/^\s+|\s+$//g; } @ele;
    die "No Masses for finding element symbol in $_\n" unless(@ele);
    @current_elements{@ele} = 1;
    my @current_elements = keys  %current_elements;
    my $species = "";
    my $starting_magnetization = "";
    my $counter = 0;
    for (@current_elements){
        $counter++;
        #my $temp = join(" ",@{$_});
        $species .= "$pot{$_}\n";
        $starting_magnetization .= "starting_magnetization($counter) = 0.01\n";
    }
    chomp ($species,$starting_magnetization);
    #my $myelement = join ('',@myelement);
    #print "$species,$starting_magnetization\n";

    my $data_path = `dirname $id`;
    my $data_name = `basename $id`;
    $data_name =~ s/\.data//g;
    chomp ($data_path, $data_name);
        open my $database ,"< $id";      
        my @data =<$database>;
        close $database;
        my %para = (
            xy => "0.0",#if no value
            xz => "0.0",
            yz => "0.0",
            natom => "",
            xlo => "",
            xhi => "",
            ylo => "",
            yhi => "",
            zlo => "",
            zhi => "",
            coords => []
        );
        for (@data){
            chomp;
    ####atoms###
            if(/(\d+)\s+atoms/){ 
                $para{natom} = $1;
            }
    ####CELL_PARAMETERS###
            elsif(/([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+xlo\s+xhi/){
                $para{xlo} = $1;
                $para{xhi} = $2;
            }
            elsif(/([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+ylo\s+yhi/){
                $para{ylo} = $1;
                $para{yhi} = $2;
            }
            elsif(/([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+zlo\s+zhi/){
                $para{zlo} = $1;
                $para{zhi} = $2;
            }
            elsif(/([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+xy\s+xz\s+yz/){
                $para{xy} = $1;
                $para{xz} = $2;
                $para{yz} = $3;
            }
#1 1 4.458517505863 1.201338326940 0.873835074284 without charge
            elsif(/\d+\s+(\d+)\s+([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)$/){
                my $ele = $ele[$1-1];
                my $x = $2 - $para{xlo};
                my $y = $3 - $para{ylo};
                my $z = $4 - $para{zlo};
                my $temp = join(" ",($ele, $x, $y, $z));
                #print "$temp No Charges\n";
                push @{$para{coords}},$temp;
            }
#1 1 1.000000 4.458517505863 1.201338326940 0.873835074284 with charge
            elsif(/\d+\s+(\d+)\s+([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)$/ && $2 != ""){
                my $ele = $ele[$1-1];
                my $x = $3 - $para{xlo};
                my $y = $4 - $para{ylo};
                my $z = $5 - $para{zlo};
                my $temp = join(" ",($ele, $x, $y, $z));
                #print "$temp have Charges\n";
                push @{$para{coords}},$temp;
            }
        }#one data file
        #check all data
        for my $k (sort keys %para){
           die "\$para{$k} is empty for $id\n" unless($para{$k});
        }
        
        my $coords = join("\n",@{$para{coords}});
        #cell parameter
        my $lx = sprintf("%.6f",$para{xhi}-$para{xlo});
        my $ly = sprintf("%.6f",$para{yhi}-$para{ylo});
        my $lz = sprintf("%.6f",$para{zhi}-$para{zlo});
        my $a = join(" ",($lx,"0.00","0.00"));
        my $b = join(" ",(sprintf("%.6f",$para{xy}),$ly,"0.00"));
        my $c = join(" ",(sprintf("%.6f",$para{xz}),sprintf("%.6f",$para{yz}),$lz));
        my $cell = join("\n",($a,$b,$c));
        my $QEoutput = $data_name .".in";
    my %QE_para = (
            calculation => "$calculation",
            output_file => "$currentPath/data2QE4MatCld/$QEoutput",
            pseudo_dir => $pseudo_dir,
            coord => $coords,
            temp => 1.0,
            press => 0.0,
            ntyp => $typeNum,
            dt => $stepsize, ###timestep size,
            nat => $para{natom},
            nstep => $nstep,
            cell_para => $cell,
            pot_spec => $species,
            starting_magnetization => $starting_magnetization,
            rho_cutoff => $rho_cutoff,
            cutoff => $cutoff
            );
      
        &QEinput(\%QE_para);
}
#all data files
  
#####here doc for QE input##########
sub QEinput
{

my ($QE_hr) = @_;

my $QEinput = <<"END_MESSAGE";
&CONTROL
calculation = "$QE_hr->{calculation}"
nstep = $QE_hr->{nstep}
etot_conv_thr = 1.0d-5
forc_conv_thr = 1.0d-4
disk_io = '/dev/null'
pseudo_dir = '$QE_hr->{pseudo_dir}'
tprnfor = .true.
tstress = .true.
verbosity = 'high'
dt = $QE_hr->{dt}
/
&SYSTEM
ntyp =  $QE_hr->{ntyp}
occupations = 'smearing' !
smearing = 'gaussian'
degauss =   0.035
ecutrho =   $rho_cutoff
ecutwfc =   $cutoff 
ibrav = 0
nat = $QE_hr->{nat}
nosym = .TRUE.
!vdw_corr = 'DFT-D3' !you may modify it using ModQEsetting.pl
$QE_hr->{starting_magnetization}
nspin = 2
/
&ELECTRONS
conv_thr =   2.d-6
electron_maxstep = 200
mixing_beta =   0.2
mixing_mode = 'plain' !'local-TF'
mixing_ndim = 8 !set 4 or 3 if OOM-killer exists (out of memory)
diagonalization = 'david' !set cg if if OOM-killer exists (out of memory). other types can be used for scf problem.
diago_david_ndim = 4 !If david is used for diagonalization. set 2 if OOM-killer appears.
/
&IONS
ion_dynamics = "beeman"
ion_temperature = "rescaling"
tempw = $QE_hr->{temp}
/
&CELL
!press_conv_thr = 0.1
cell_dynamics = "pr"
press = $QE_hr->{press}
cell_dofree = "all"
/
K_POINTS {automatic}
2 2 2 0 0 0
ATOMIC_SPECIES
$QE_hr->{pot_spec}
ATOMIC_POSITIONS {angstrom}
$QE_hr->{coord}
CELL_PARAMETERS {angstrom}
$QE_hr->{cell_para}
END_MESSAGE

my $file = $QE_hr->{output_file};
open(FH, '>', $QE_hr->{output_file}) or die $!;
print FH $QEinput;
close(FH);
}
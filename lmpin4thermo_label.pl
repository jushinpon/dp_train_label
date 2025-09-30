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

###parameters to set first
my $currentPath = getcwd();# dir for all scripts
chdir("..");
my $mainPath = getcwd();# main path of Perl4dpgen dir
chdir("$currentPath");

my @datafile = `find $currentPath/data4thermo -name "*.data"|grep -v 2d`;#find all data files
map { s/^\s+|\s+$//g; } @datafile;

my $temperatur_initial = 300;
my $temperatur_end = 1000;
my @pressure = ("0");#10000 in lammps for 1 Gpa
my $run_step = 500000;
my $timestep = 0.001;
my $tdamp = $timestep*100;
my $pdamp = $timestep*1000;
my $out_freq = 500;#dlp md.out frequency
###minimize set

my $PotentialPath = "$mainPath/dp_train_new/dp_train";
my @pb_files = `find $PotentialPath -type f -name "*.pb"|grep compress|grep -v -- "-p"`;#all npy files
map { s/^\s+|\s+$//g; } @pb_files;
die "No DLP pb files \n" unless (@pb_files);
my $Potential_pb = join (" ",@pb_files);

my $Potential_prod="$Potential_pb out_file md.out out_freq $out_freq";
#for surface $box_relax,should be modified
my $box_relax = "z 0.0";#"aniso 200.0";#"x 200 y 200 z 200";#"aniso 0.0"; # a little compressed to get some compressed structures
my $min_value = "0.0 0.0 50000 100000"; #etol ftol maxiter maxeval
###NVT (0) or NPT (1) set
my $ensemble = "1" ; # If the ensemble is npt or nvt. Set according to your own needs.
my $ensemble_name;

if ($ensemble == 0){
    $ensemble_name = "NVT";
} else {
    $ensemble_name = "NPT";
} # # If the $ensemble is npt, change $ensemble_name to NPT. If the $ensemble is nvt, $ensemble_name is changed to NVT.
## for here doc
`rm -rf $currentPath/thermo_label`;
`mkdir $currentPath/thermo_label`;

for my $id (@datafile){
    my $data_path = `dirname $id`;
    $data_path =~ s/^\s+|\s+$//g;
    my $data_name = `basename $id`;
    $data_name =~ s/^\s+|\s+$//g;

    $data_name =~ s/\.data//g;
    for my $press (@pressure){
        my $lmpoutput = "$data_name-lmpT$temperatur_initial"."to$temperatur_end"."-P$press"."-R$run_step.in";
        my $file_dir = "$currentPath/thermo_label/$data_name-lmpT$temperatur_initial"."to$temperatur_end"."-P$press"."-R$run_step";
        `mkdir -p $file_dir`;
        #print "file_dir: $file_dir\n";
        
        my %lmp_para = (
            read_data => $id,
            Potential_min => "$pb_files[0]",#use the first one
            Potential_prod => "$Potential_prod",# for label
            output_file => "$file_dir/$lmpoutput",
            infile =>$lmpoutput,
            press => $press,
            ts => $timestep,
            tdamp => $tdamp,
            pdamp => $pdamp,
            min_value => $min_value,
            temperatur_initial => $temperatur_initial,
            temperatur_end => $temperatur_end,
            ensemble => $ensemble,
            ensemble_name =>$ensemble_name,
            run_step => $run_step,
            out_freq => $out_freq,
            currentPath => $currentPath,
            data_name => $data_name,
            box_relax => $box_relax,
            );
      
        &lmpinput(\%lmp_para);
    }#pressure
    
}#all data files
   
######here doc for LAMMPS input##########
sub lmpinput
{

my ($lmp_hr) = @_;

my $lmpinput = <<"END_MESSAGE";
#echo none
units metal 
dimension 3 
boundary p p p 
box tilt large
atom_style atomic 
atom_modify map array 
shell mkdir lmp_output
# ---------- Create Atoms --------------------- 
read_data $lmp_hr->{read_data}
# ---------- Define Interatomic Potential --------------------- 
pair_style deepmd $lmp_hr->{Potential_min}
pair_coeff * *  
#----------------------------------------------
neighbor 1 bin 
neigh_modify delay 10 every 5 check yes one 5000
#-----------------minimize---------------------
min_style cg
fix 5 all box/relax $lmp_hr->{box_relax}
timestep $lmp_hr->{ts}
min_style	     cg
thermo 100
thermo_style custom step temp density lx press pxx pyy pzz pxy pxz pyz pe
#dump 1 all custom $lmp_hr->{out_freq} MIN_*.cfg id type x y z xu yu zu
minimize $lmp_hr->{min_value}
unfix 5
#undump 1
write_data $lmp_hr->{data_name}_minimize.data
clear
#------------------ensemble-------------------------
units metal 
dimension 3 
boundary p p p 
box tilt large
atom_style atomic 
atom_modify map array
variable ensemble equal $lmp_hr->{ensemble}
# ---------- Create Atoms --------------------- 
read_data $lmp_hr->{data_name}_minimize.data
# ---------- Define Interatomic Potential --------------------- 
pair_style deepmd $lmp_hr->{Potential_prod}
pair_coeff * *
#----------------------------------------------
neighbor 1 bin 
neigh_modify delay 10 every 5 check yes one 5000
#----------------------------------------------
reset_timestep 0
timestep $lmp_hr->{ts}
velocity all create $lmp_hr->{temperatur_initial} 12345 mom yes rot yes dist gaussian
velocity all scale $lmp_hr->{temperatur_initial}
if " \${ensemble} == 0 " then &
"fix 1 all nvt temp $lmp_hr->{temperatur_initial} $lmp_hr->{temperatur_end} $lmp_hr->{tdamp}" &
else &
"fix 1 all npt temp $lmp_hr->{temperatur_initial} $lmp_hr->{temperatur_end} $lmp_hr->{tdamp} z $lmp_hr->{press} $lmp_hr->{press} $lmp_hr->{pdamp}"
#"fix 1 all npt temp $lmp_hr->{temperatur_initial} $lmp_hr->{temperatur_end} $lmp_hr->{tdamp} aniso $lmp_hr->{press} $lmp_hr->{press} $lmp_hr->{pdamp}"

thermo 100
thermo_style custom step temp density pxx pyy pzz press pe
shell cd lmp_output
dump 1 all custom $lmp_hr->{out_freq} lmp_*.cfg id type x y z xu yu zu
run $lmp_hr->{run_step}
unfix 1
undump 1
shell cd ..
write_data $lmp_hr->{data_name}_$lmp_hr->{ensemble_name}.data

END_MESSAGE

my $file = $lmp_hr->{output_file};
#print "myfile: $file\n";
open(FH, '>', "$file") or die $!;
print FH $lmpinput;
close(FH);
}
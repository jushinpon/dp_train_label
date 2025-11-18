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

my $PotentialPath = "$mainPath/dp_train_new/dp_train";

my @datafile = `find $currentPath/data4shear -name "*.data"`;#find all data files
map { s/^\s+|\s+$//g; } @datafile;

#my $temperatur_initial = 300;
#my $temperatur_end = 1000;
my @pressure = ("0");#10000 in lammps for 1 Gpa
#my $run_step = 500000;
#my $timestep = "0.001";
#my $tdamp = "0.1";
#my $pdamp = "1.0";

my $out_freq = 1;#dlp md.out frequency
my $loop = 21; #total loop steps

# shear range: gamma from -0.02 to 0.02
my $min_scale = -0.03;#for shear strain min (gamma)
my $max_scale = 0.03;#for shear strain max (gamma)

my @pb_files = `find $PotentialPath -type f -name "*.pb"|grep compress|grep -v -- "-p"`;#all npy files
map { s/^\s+|\s+$//g; } @pb_files;

my $Potential_pb = join (" ",@pb_files);
my $Potential_prod="$Potential_pb out_file md.out out_freq $out_freq";

#for surface $box_relax,should be modified
my $box_relax = "iso 0.0 nreset 50";#"aniso 200.0";#"x 200 y 200...so 0.0"; # a little compressed to get some compressed structures
my $min_value = "0.0 0.0 50000 100000"; #etol ftol maxiter maxeval

#my $ensemble = "1" ; # If the ensemble is npt or nvt. Set according to your own needs.
#my $ensemble_name;
#
#if ($ensemble == 0){
#    $ensemble_name = "NVT";
#} else {
#    $ensemble_name = "NPT";
#} # # If the $ensemble is npt, change $ensemble_name to NPT. If the $ensemble is nvt, $ensemble_name is changed to NVT.
## for here doc
`rm -rf $currentPath/shear_label`;
`mkdir $currentPath/shear_label`;

for my $id (@datafile){
    my $data_path = `dirname $id`;
    $data_path =~ s/^\s+|\s+$//g;
    my $data_name = `basename $id`;
    $data_name =~ s/^\s+|\s+$//g;

    $data_name =~ s/\.data//g;
    for my $press (@pressure){
        my $lmpoutput = "$data_name-lmpshear.in";
        my $file_dir = "$currentPath/shear_label/$data_name-lmpshear";
        `mkdir -p $file_dir`;

        my %lmp_para = (
            read_data => $id,
            Potential_min => "$pb_files[0]",#use the first one
            Potential_prod => "$Potential_prod",# for label
            output_file => "$file_dir/$lmpoutput",
            infile =>$lmpoutput,
            #press => $press,
            #ts => $timestep,
            #tdamp => $tdamp,
            #pdamp => $pdamp,
            min_value => $min_value,
            #temperatur_initial => $temperatur_initial,
            #temperatur_end => $temperatur_end,
            #ensemble => $ensemble,
            #ensemble_name =>$ensemble_name,
            #run_step => $run_step,
            loop => $loop,
            min_scale => $min_scale,
            max_scale => $max_scale,
            out_freq => $out_freq,
            currentPath => $currentPath,
            data_name => $data_name,
            box_relax => $box_relax,
            );
  
        &lmpinput(\%lmp_para);
    }
}

sub lmpinput
{

my ($lmp_hr) = @_;

my $lmpinput = <<"END_MESSAGE";
#echo none
plugin load libdeepmd_lmp.so #for deepmd-v3 only
units metal 
dimension 3 
boundary p p p 
box tilt large
atom_style atomic 
atom_modify map array 
shell mkdir lmp_output

# ---------- Create Atoms (equilibrium minimize) ---------------------
read_data $lmp_hr->{read_data}

# ---------- Define Interatomic Potential (minimize) -----------------
pair_style deepmd $lmp_hr->{Potential_min}
pair_coeff * *  
#----------------------------------------------
neighbor 1 bin 
neigh_modify delay 10 every 5 check yes one 5000
#-----------------minimize---------------------
min_style cg
fix freeze all setforce 0.0 0.0 0.0
fix 5 all box/relax $lmp_hr->{box_relax}
#timestep $lmp_hr->{ts}
min_style	     cg
thermo 100
thermo_style custom step temp density lx press pxx pyy pzz pxy pxz pyz pe
dump 1 all custom 50 MIN_*.cfg id type x y z xu yu zu
#dump 1 all custom $lmp_hr->{out_freq} MIN_*.cfg id type x y z xu yu zu
minimize $lmp_hr->{min_value}
unfix 5
unfix freeze
undump 1
write_data $lmp_hr->{data_name}_minimize.data

clear
#------------------shear ensemble-------------------------
plugin load libdeepmd_lmp.so #for deepmd-v3 only
units metal 
dimension 3 
boundary p p p 
box tilt large
atom_style atomic 
atom_modify map array
# ---------- Create Atoms (start from minimized structure) ----------
read_data $lmp_hr->{data_name}_minimize.data

# --- store initial box lengths (snapshot) ---
variable lx equal lx
variable ly equal ly
variable lz equal lz
variable lx0 equal \${lx}
variable ly0 equal \${ly}
variable lz0 equal \${lz}
variable lx0 equal \${lx0}
variable ly0 equal \${ly0}
variable lz0 equal \${lz0}

# --- store initial tilt factors (snapshot, usually 0) ---
variable xy0 equal xy
variable xz0 equal xz
variable yz0 equal yz
variable xy0 equal \${xy0}
variable xz0 equal \${xz0}
variable yz0 equal \${yz0}

# --- number of steps per shear series (e.g. 21) ---
variable nstep equal $lmp_hr->{loop}

# --- shear range: min_scale/max_scale are +/- gamma ---
# gamma_xy = xy / Ly  ->  xy = gamma_xy * Ly
# gamma_xz = xz / Lz  ->  xz = gamma_xz * Lz
# gamma_yz = yz / Lz  ->  yz = yz / Lz
variable xy_min equal \${ly0}*$lmp_hr->{min_scale}
variable xy_max equal \${ly0}*$lmp_hr->{max_scale}
variable xz_min equal \${lz0}*$lmp_hr->{min_scale}
variable xz_max equal \${lz0}*$lmp_hr->{max_scale}
variable yz_min equal \${lz0}*$lmp_hr->{min_scale}
variable yz_max equal \${lz0}*$lmp_hr->{max_scale}

variable dxy equal (\${xy_max}-\${xy_min})/(v_nstep-1)
variable dxz equal (\${xz_max}-\${xz_min})/(v_nstep-1)
variable dyz equal (\${yz_max}-\${yz_min})/(v_nstep-1)

# ---------- Define Interatomic Potential (production/shear) -------
pair_style deepmd $lmp_hr->{Potential_prod}
pair_coeff * * 
#----------------------------------------------
neighbor 1 bin
neigh_modify delay 1 every 1 check yes one 5000
#----------------------------------------------
reset_timestep 0

thermo 1
thermo_style custom step temp density pxx pyy pzz pxy pxz pyz press pe

# ============================================================
# 1) xy shear: nstep steps, block = 0
# ============================================================

# reset to unreformed box before xy shear
change_box all x final 0.0 \${lx0} y final 0.0 \${ly0} z final 0.0 \${lz0} xy final \${xy0} xz final \${xz0} yz final \${yz0} remap units box

variable i loop \${nstep}
label loop_xy

  variable ci equal \${i}-1      # index: 0 .. nstep-1
  variable block equal 0
  variable step_index equal "v_block*v_nstep + v_ci"
  reset_timestep \${step_index}

  variable new_xy equal "v_xy_min + v_dxy*v_ci"

  change_box all x  final 0.0 \${lx0} y  final 0.0 \${ly0} z  final 0.0 \${lz0} xy final \${new_xy} xz final \${xz0} yz final \${yz0} remap units box

  shell cd lmp_output
  dump 1 all custom $lmp_hr->{out_freq} lmp_*.cfg id type x y z xu yu zu
  run 0
  undump 1
  shell cd ..

  next i
jump SELF loop_xy

# ============================================================
# 2) xz shear: nstep steps, block = 1
# ============================================================

# reset to unreformed box before xz shear
change_box all x final 0.0 \${lx0} y final 0.0 \${ly0} z final 0.0 \${lz0} xy final \${xy0} xz final \${xz0} yz final \${yz0} remap units box

variable j loop \${nstep}
label loop_xz

  variable cj equal \${j}-1      # index: 0 .. nstep-1
  variable block equal 1
  variable step_index equal "v_block*v_nstep + v_cj"
  reset_timestep \${step_index}

  variable new_xz equal "v_xz_min + v_dxz*v_cj"

  change_box all x  final 0.0 \${lx0} y  final 0.0 \${ly0} z  final 0.0 \${lz0} xy final \${xy0} xz final \${new_xz} yz final \${yz0} remap units box

  shell cd lmp_output
  dump 1 all custom $lmp_hr->{out_freq} lmp_*.cfg id type x y z xu yu zu
  run 0
  undump 1
  shell cd ..

  next j
jump SELF loop_xz

# ============================================================
# 3) yz shear: nstep steps, block = 2
# ============================================================

# reset to unreformed box before yz shear
change_box all x final 0.0 \${lx0} y final 0.0 \${ly0} z final 0.0 \${lz0} xy final \${xy0} xz final \${xz0} yz final \${yz0} remap units box

variable k loop \${nstep}
label loop_yz

  variable ck equal \${k}-1      # index: 0 .. nstep-1
  variable block equal 2
  variable step_index equal "v_block*v_nstep + v_ck"
  reset_timestep \${step_index}

  variable new_yz equal "v_yz_min + v_dyz*v_ck"

  change_box all x  final 0.0 \${lx0} y  final 0.0 \${ly0} z  final 0.0 \${lz0} xy final \${xy0} xz final \${xz0} yz final \${new_yz} remap units box

  shell cd lmp_output
  dump 1 all custom $lmp_hr->{out_freq} lmp_*.cfg id type x y z xu yu zu
  run 0
  undump 1
  shell cd ..

  next k
jump SELF loop_yz

write_data $lmp_hr->{data_name}_final.data

END_MESSAGE

my $file = $lmp_hr->{output_file};
#print "myfile: $file\n";
open(FH, '>', "$file") or die $!;
print FH $lmpinput;
close(FH);
}

#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(ceil floor);
use File::Basename;
use File::Path qw(make_path);

# ==============================================================================
# USER CONFIGURATION
# ==============================================================================
my $data_sour      = "/home/jsp/SnPbTe_alloys/dp_train_new/initial/Sn15Pb17Te32-T300-P0/Sn15Pb17Te32-T300-P0.data";
my $QEinput_source = "/home/jsp/SnPbTe_alloys/dp_train_new/initial/Sn15Pb17Te32-T300-P0/Sn15Pb17Te32-T300-P0.in";
my @xyz_supercell  = (2, 2, 2); # supercell in x, y, and z dimensions

# ==============================================================================
# 1. SETUP OUTPUT DIRECTORY
# ==============================================================================
my $basename = basename($QEinput_source, ".in");
my $super_str = join("", @xyz_supercell); 
my $out_dir = "super${super_str}-${basename}";

unless (-d $out_dir) {
    make_path($out_dir) or die "Failed to create path: $out_dir";
    print "Created directory: $out_dir\n";
}

my $output_qe   = "$out_dir/${basename}.in";
my $output_data = "$out_dir/${basename}.data";

# ==============================================================================
# 2. PARSE QE INPUT (Species & K-Points)
# ==============================================================================
open(my $fh_in, '<', $QEinput_source) or die "Could not open QE input: $!";
my @qe_lines = <$fh_in>;
close($fh_in);

my @species_map = (); 
my @original_k = ();
my $found_k = 0;
my $k_type = "";
my $in_species = 0;

foreach my $line (@qe_lines) {
    if ($line =~ /ATOMIC_SPECIES/i) { $in_species = 1; next; }
    
    if ($in_species) {
        if ($line =~ /^\s*&(CONTROL|SYSTEM|ELECTRONS|IONS|CELL)/i || 
            $line =~ /ATOMIC_POSITIONS/i || 
            $line =~ /K_POINTS/i || 
            $line =~ /CELL_PARAMETERS/i) {
            $in_species = 0;
        }
    }
    
    if ($in_species && $line =~ /\S/) {
        my ($label) = split(/\s+/, $line =~ s/^\s+//r);
        push(@species_map, $label) if ($label && $label !~ /^!/);
    }
}

for (my $i = 0; $i < scalar(@qe_lines); $i++) {
    if ($qe_lines[$i] =~ /K_POINTS\s*\{?(\w+)\}?/i) {
        $k_type = lc($1);
        if ($k_type eq "automatic") {
            my $k_line = $qe_lines[$i+1];
            $k_line =~ s/^\s+|\s+$//g;
            @original_k = split(/\s+/, $k_line);
            $found_k = 1;
        } else {
            die "Error: Script only supports 'K_POINTS automatic'. Found: $k_type\n";
        }
        last;
    }
}

if (!@species_map) { die "Error: Could not find ATOMIC_SPECIES.\n"; }
if (!$found_k) { die "Error: Could not find K_POINTS automatic.\n"; }

# ==============================================================================
# 3. PARSE LAMMPS DATA (Geometry & Masses)
# ==============================================================================
open(my $fh_data, '<', $data_sour) or die "Could not open Data file: $!";
my @data_lines = <$fh_data>;
close($fh_data);

my ($xlo, $xhi, $ylo, $yhi, $zlo, $zhi) = (0,0,0,0,0,0);
my ($xy, $xz, $yz) = (0,0,0);
my @atoms = ();
my @masses = ();
my $in_atoms = 0;
my $in_masses = 0;

foreach my $line (@data_lines) {
    # Parse Box Info
    if ($line =~ /xlo\s+xhi/) { ($xlo, $xhi) = split(/\s+/, $line); }
    if ($line =~ /ylo\s+yhi/) { ($ylo, $yhi) = split(/\s+/, $line); }
    if ($line =~ /zlo\s+zhi/) { ($zlo, $zhi) = split(/\s+/, $line); }
    if ($line =~ /xy\s+xz\s+yz/) { ($xy, $xz, $yz) = split(/\s+/, $line); }

    # Parse Masses (Store them to write back later)
    if ($line =~ /^Masses/) { $in_masses = 1; $in_atoms = 0; next; }
    if ($in_masses && $line =~ /^\s*[A-Za-z]/) { $in_masses = 0; } # Stop if next section name
    if ($in_masses && $line =~ /^\s*\d+/) { push(@masses, $line); }

    # Parse Atoms
    if ($line =~ /^Atoms/) { $in_atoms = 1; $in_masses = 0; next; }
    if ($in_atoms && $line =~ /^\s*$/) { next; }
    if ($in_atoms && $line =~ /^\s*\d+/) {
        my @parts = split(/\s+/, $line =~ s/^\s+//r);
        my $z = pop(@parts);
        my $y = pop(@parts);
        my $x = pop(@parts);
        my $type = $parts[1]; 
        push(@atoms, { type => $type, x => $x, y => $y, z => $z });
    }
}

# Basis vectors (Original)
my $lx = $xhi - $xlo;
my $ly = $yhi - $ylo;
my $lz = $zhi - $zlo;
my @vec_a = ($lx, 0, 0);
my @vec_b = ($xy, $ly, 0);
my @vec_c = ($xz, $yz, $lz);

# ==============================================================================
# 4. GENERATE SUPERCELL
# ==============================================================================
my ($sx, $sy, $sz) = @xyz_supercell;
my @new_atoms = ();

for my $ix (0 .. $sx - 1) {
    for my $iy (0 .. $sy - 1) {
        for my $iz (0 .. $sz - 1) {
            my $shift_x = $ix * $vec_a[0] + $iy * $vec_b[0] + $iz * $vec_c[0];
            my $shift_y = $ix * $vec_a[1] + $iy * $vec_b[1] + $iz * $vec_c[1];
            my $shift_z = $ix * $vec_a[2] + $iy * $vec_b[2] + $iz * $vec_c[2];

            foreach my $atom (@atoms) {
                push(@new_atoms, {
                    type => $atom->{type},
                    x    => $atom->{x} + $shift_x,
                    y    => $atom->{y} + $shift_y,
                    z    => $atom->{z} + $shift_z
                });
            }
        }
    }
}

# Scale Lattice Vectors for QE
my @new_vec_a = ($vec_a[0]*$sx, $vec_a[1]*$sx, $vec_a[2]*$sx);
my @new_vec_b = ($vec_b[0]*$sy, $vec_b[1]*$sy, $vec_b[2]*$sy);
my @new_vec_c = ($vec_c[0]*$sz, $vec_c[1]*$sz, $vec_c[2]*$sz);

# Scale K-Points
my @new_k = @original_k;
$new_k[0] = ceil($original_k[0] / $sx);
$new_k[1] = ceil($original_k[1] / $sy);
$new_k[2] = ceil($original_k[2] / $sz);
$new_k[0] = 1 if $new_k[0] < 1;
$new_k[1] = 1 if $new_k[1] < 1;
$new_k[2] = 1 if $new_k[2] < 1;

# ==============================================================================
# 5. WRITE QE OUTPUT (.in)
# ==============================================================================
open(my $fh_out, '>', $output_qe) or die "Could not write QE output: $!";
print "Writing QE input to: $output_qe\n";

my $skip_mode = 0;

foreach my $line (@qe_lines) {
    # Fix NAT format: Start of line, no comma
    if ($line =~ /nat\s*=/) {
        print $fh_out "nat = " . scalar(@new_atoms) . "\n";
        next;
    }

    if ($line =~ /ATOMIC_POSITIONS/i) {
        $skip_mode = 1;
        print $fh_out "ATOMIC_POSITIONS {angstrom}\n";
        foreach my $atom (@new_atoms) {
            my $label = $species_map[$atom->{type} - 1] // "X";
            printf $fh_out "%-4s  %16.10f  %16.10f  %16.10f\n", $label, $atom->{x}, $atom->{y}, $atom->{z};
        }
        next;
    }
    
    if ($line =~ /CELL_PARAMETERS/i) {
        $skip_mode = 1;
        print $fh_out "CELL_PARAMETERS {angstrom}\n";
        printf $fh_out "  %16.10f  %16.10f  %16.10f\n", @new_vec_a;
        printf $fh_out "  %16.10f  %16.10f  %16.10f\n", @new_vec_b;
        printf $fh_out "  %16.10f  %16.10f  %16.10f\n", @new_vec_c;
        next;
    }

    if ($line =~ /K_POINTS/i) {
        $skip_mode = 1;
        print $fh_out "K_POINTS {automatic}\n";
        printf $fh_out "  %d %d %d %d %d %d\n", @new_k;
        next;
    }

    if ($skip_mode) {
        if ($line =~ /^\s*&(CONTROL|SYSTEM|ELECTRONS|IONS|CELL)/i || 
            $line =~ /ATOMIC_SPECIES/i ||
            $line =~ /ATOMIC_POSITIONS/i || 
            $line =~ /K_POINTS/i || 
            $line =~ /CELL_PARAMETERS/i) {
             # If line is not just data (numbers), stop skipping
             unless ($line =~ /^\s*[0-9.-]+\s+[0-9.-]+\s+[0-9.-]+/ && $line !~ /[a-zA-Z]/) {
                 $skip_mode = 0;
             }
        }
    }
    print $fh_out $line unless $skip_mode;
}
close($fh_out);

# ==============================================================================
# 6. WRITE LAMMPS DATA OUTPUT (.data)
# ==============================================================================
open(my $fh_dat, '>', $output_data) or die "Could not write Data output: $!";
print "Writing LAMMPS data to: $output_data\n";

# New Box Dimensions
my $lx_new = $lx * $sx;
my $ly_new = $ly * $sy;
my $lz_new = $lz * $sz;
my $xy_new = $xy * $sy;
my $xz_new = $xz * $sz;
my $yz_new = $yz * $sz;

print $fh_dat "LAMMPS data file generated by script. Supercell $sx x $sy x $sz\n\n";
print $fh_dat scalar(@new_atoms) . " atoms\n";
print $fh_dat scalar(@species_map) . " atom types\n\n";

printf $fh_dat "0.0 %16.10f xlo xhi\n", $lx_new;
printf $fh_dat "0.0 %16.10f ylo yhi\n", $ly_new;
printf $fh_dat "0.0 %16.10f zlo zhi\n", $lz_new;

if ($xy != 0 || $xz != 0 || $yz != 0) {
    printf $fh_dat "%16.10f %16.10f %16.10f xy xz yz\n", $xy_new, $xz_new, $yz_new;
}
print $fh_dat "\n";

if (@masses) {
    print $fh_dat "Masses\n\n";
    print $fh_dat @masses;
    print $fh_dat "\n";
}

print $fh_dat "Atoms # atomic\n\n";
my $atom_id = 1;
foreach my $atom (@new_atoms) {
    # ID type x y z
    printf $fh_dat "%d %d %16.10f %16.10f %16.10f\n", $atom_id, $atom->{type}, $atom->{x}, $atom->{y}, $atom->{z};
    $atom_id++;
}
close($fh_dat);

print "Processing Complete.\n";
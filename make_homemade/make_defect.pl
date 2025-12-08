#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(shuffle);
use File::Basename;
use Cwd 'abs_path';

# --- User Settings ---
# Please provide the absolute path to your source files here
my $data_sour = "/home/jsp/SnPbTe_alloys/dp_train_new/initial/Sn15Pb17Te32-T300-P0/Sn15Pb17Te32-T300-P0.data";
my $QEinput_source = "/home/jsp/SnPbTe_alloys/dp_train_new/initial/Sn15Pb17Te32-T300-P0/Sn15Pb17Te32-T300-P0.in";

my @delete_element = ("Sn", "Pb", "Te", "Te"); # Element type to be deleted
my @delete_num     = (1, 1, 1, 2);              # Corresponding number to be deleted

# --- Main Processing ---

# Resolve absolute paths to ensure we can find files from anywhere
# This handles cases where you might use relative paths like ../file.data
my $abs_data_sour = abs_path($data_sour);
my $abs_qe_source = abs_path($QEinput_source);

unless (-e $abs_data_sour && -e $abs_qe_source) {
    die "Error: One or more input files not found:\nData: $data_sour\nQE: $QEinput_source\n";
}

# Parse base filename for naming the new folders
# e.g., /path/to/Sn15Pb17Te32-T300-P0.data -> Sn15Pb17Te32-T300-P0
my ($base_name, $dir, $suffix) = fileparse($abs_data_sour, qr/\.[^.]*/);

print "Processing base case: $base_name\n";

# 1. Read and Parse LAMMPS Data File
my %lammps_info = read_lammps_data($abs_data_sour);

# 2. Read and Parse QE Input File
my %qe_info = read_qe_input($abs_qe_source);

# Check consistency
if ($lammps_info{total_atoms} != $qe_info{total_atoms}) {
    die "Error: Atom count mismatch! LAMMPS ($lammps_info{total_atoms}) vs QE ($qe_info{total_atoms}).\n";
}

# 3. Loop through the defect tasks
for (my $i = 0; $i < scalar(@delete_element); $i++) {
    my $elem = $delete_element[$i];
    my $num  = $delete_num[$i];
    
    # Construct folder name prefix (e.g., def-Sn1)
    my $prefix = "def-${elem}${num}";
    my $folder_name = "${prefix}-${base_name}";
    
    # Define absolute path for the new folder (created in the current working directory)
    my $current_dir = `pwd`;
    chomp($current_dir);
    my $target_dir = "$current_dir/$folder_name";
    
    print "\n--- Case $i: Creating $target_dir ---\n";

    # --- Use Shell Command to Create Directory ---
    if (-d $target_dir) {
        `rm -rf "$target_dir"`; # Clean existing folder
    }
    `mkdir -p "$target_dir"`;

    # Identify Atom Type ID for the element (e.g., Sn -> 1)
    my $type_id = $lammps_info{element_map}{$elem};
    unless (defined $type_id) {
        warn "Warning: Element '$elem' not found in LAMMPS Masses. Skipping.\n";
        next;
    }

    # Find candidates
    my @candidate_indices;
    for (my $idx = 0; $idx < $lammps_info{total_atoms}; $idx++) {
        # lammps_atoms: [id, type, x, y, z, line_content]
        if ($lammps_info{atoms}[$idx]->[1] == $type_id) {
            push @candidate_indices, $idx;
        }
    }

    if (scalar(@candidate_indices) < $num) {
        warn "Not enough atoms to delete. Skipping.\n";
        next;
    }

    # Randomly select atoms to delete
    my @shuffled = shuffle(@candidate_indices);
    my @to_delete = @shuffled[0 .. $num-1];
    
    # Hash for fast lookup of deleted indices
    my %deleted_lookup = map { $_ => 1 } @to_delete;

    # --- Write New Files ---
    my $out_data = "$target_dir/${folder_name}.data";
    my $out_qe   = "$target_dir/${folder_name}.in";
    
    # 1. Write LAMMPS Data
    open(my $fh_d, '>', $out_data) or die "Cannot write $out_data: $!";
    
    # Adjust header atom count
    my $new_total = $lammps_info{total_atoms} - $num;
    foreach my $line (@{$lammps_info{header}}) {
        if ($line =~ /(\d+)\s+atoms/) {
            print $fh_d "$new_total atoms\n";
        } else {
            print $fh_d $line;
        }
    }
    
    # Write atoms with renumbered IDs
    my $new_id = 1;
    for (my $idx = 0; $idx < $lammps_info{total_atoms}; $idx++) {
        next if $deleted_lookup{$idx}; # Skip deleted

        my $line = $lammps_info{atoms}[$idx]->[5];
        # Replace the old ID (first number) with the new continuous ID
        $line =~ s/^\s*(\d+)/$new_id/;
        print $fh_d $line;
        $new_id++;
    }
    close $fh_d;

    # 2. Write QE Input
    open(my $fh_q, '>', $out_qe) or die "Cannot write $out_qe: $!";
    
    # Write Pre-Atoms part (update nat)
    foreach my $line (@{$qe_info{pre_atoms}}) {
        if ($line =~ /nat\s*=\s*\d+/) {
            $line =~ s/nat\s*=\s*\d+/nat = $new_total/;
        }
        print $fh_q $line;
    }
    
    # Write Atoms (skip deleted)
    for (my $idx = 0; $idx < $qe_info{total_atoms}; $idx++) {
        next if $deleted_lookup{$idx};
        print $fh_q $qe_info{atom_lines}[$idx];
    }
    
    # Write Post-Atoms part
    print $fh_q $_ for @{$qe_info{post_atoms}};
    close $fh_q;

    print "Created: $out_data\n";
    print "Created: $out_qe\n";
}

print "\nDone.\n";

# --- Subroutines for Parsing ---

sub read_lammps_data {
    my ($file) = @_;
    open(my $fh, '<', $file) or die "Cannot open $file: $!";
    my @lines = <$fh>;
    close $fh;

    my %data;
    my (@header, @atoms);
    my %elem_map;
    my $section = "header";
    
    foreach my $line (@lines) {
        if ($line =~ /(\d+)\s+atoms/) { $data{total_atoms} = $1; }

        if ($line =~ /^Masses/) {
            $section = "masses";
            push @header, $line; next;
        } elsif ($line =~ /^Atoms/) {
            $section = "atoms_pre";
            push @header, $line; next;
        }

        if ($section eq "header") {
            push @header, $line;
        } elsif ($section eq "masses") {
            push @header, $line;
            # Parse "1 118.71 # Sn" to map Sn -> 1
            if ($line =~ /^\s*(\d+).+#\s*(\w+)/) { $elem_map{$2} = $1; }
        } elsif ($section eq "atoms_pre") {
            push @header, $line;
            # Detect start of atom data
            if ($line =~ /^\s*\d+/) { $section = "atoms_body"; redo; }
        } elsif ($section eq "atoms_body") {
            # ID Type X Y Z
            if ($line =~ /^\s*(\d+)\s+(\d+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)/) {
                push @atoms, [$1, $2, $3, $4, $5, $line];
            }
        }
    }
    $data{header} = \@header;
    $data{atoms} = \@atoms;
    $data{element_map} = \%elem_map;
    return %data;
}

sub read_qe_input {
    my ($file) = @_;
    open(my $fh, '<', $file) or die "Cannot open $file: $!";
    my @lines = <$fh>;
    close $fh;

    my %data;
    my (@pre, @atoms, @post);
    my $section = "pre";

    foreach my $line (@lines) {
        if ($section eq "pre") {
            push @pre, $line;
            if ($line =~ /ATOMIC_POSITIONS/) { $section = "atoms"; }
        } elsif ($section eq "atoms") {
            # End of atoms block if we see new keywords or slashes
            if ($line =~ /K_POINTS|CELL_PARAMETERS|ATOMIC_SPECIES|^\s*[\&\/]/) {
                $section = "post";
                push @post, $line;
            } else {
                if ($line =~ /\w/) { push @atoms, $line; } # Skip empty lines
            }
        } else {
            push @post, $line;
        }
    }
    $data{pre_atoms} = \@pre;
    $data{atom_lines} = \@atoms;
    $data{post_atoms} = \@post;
    $data{total_atoms} = scalar(@atoms);
    return %data;
}
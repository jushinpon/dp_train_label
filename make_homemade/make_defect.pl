#!/usr/bin/perl
#rm -rf /home/jsp/SnPbTe_alloys/dp_train_new/initial/def-*   ---> be very careful when using this command
#cp -r def-* /home/jsp/SnPbTe_alloys/dp_train_new/initial/

use strict;
use warnings;
use List::Util qw(shuffle);
use File::Basename;
use Cwd 'abs_path';

# --- User Settings ---
my $data_sour = "/home/jsp/SnPbTe_alloys/dp_train_new/initial/Sn15Pb17Te32-T300-P0/Sn15Pb17Te32-T300-P0.data";
my $QEinput_source = "/home/jsp/SnPbTe_alloys/dp_train_new/initial/Sn15Pb17Te32-T300-P0/Sn15Pb17Te32-T300-P0.in";

my @delete_element = ("Sn", "Pb", "Te", "Te"); # Element type to be deleted
my @delete_num     = (1, 1, 1, 2);              # Corresponding number to be deleted

# --- Main Processing ---

my $abs_data_sour = abs_path($data_sour);
my $abs_qe_source = abs_path($QEinput_source);

unless (-e $abs_data_sour && -e $abs_qe_source) {
    die "Error: One or more input files not found:\nData: $data_sour\nQE: $QEinput_source\n";
}

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
    
    my $prefix = "def-${elem}${num}";
    my $folder_name = "${prefix}-${base_name}";
    
    my $current_dir = `pwd`;
    chomp($current_dir);
    my $target_dir = "$current_dir/$folder_name";
    
    print "\n======================================================\n";
    print "Task $i: Deleting $num atom(s) of element '$elem'\n";
    print "Output Directory: $target_dir\n";

    if (-d $target_dir) {
        `rm -rf "$target_dir"`; 
    }
    `mkdir -p "$target_dir"`;

    my $type_id = $lammps_info{element_map}{$elem};
    unless (defined $type_id) {
        warn "Warning: Element '$elem' not found in LAMMPS Masses. Skipping.\n";
        next;
    }

    my @candidate_indices;
    for (my $idx = 0; $idx < $lammps_info{total_atoms}; $idx++) {
        if ($lammps_info{atoms}[$idx]->[1] == $type_id) {
            push @candidate_indices, $idx;
        }
    }

    if (scalar(@candidate_indices) < $num) {
        warn "Not enough atoms to delete. Skipping.\n";
        next;
    }

    my @shuffled = shuffle(@candidate_indices);
    my @to_delete = @shuffled[0 .. $num-1];
    
    # --- PRINT DELETED ATOM INFO ---
    print "------------------------------------------------------\n";
    print "Deleted Atoms Information:\n";
    foreach my $idx (@to_delete) {
        my $l_id   = $lammps_info{atoms}[$idx]->[0];
        my $l_type = $lammps_info{atoms}[$idx]->[1];
        my $q_line = $qe_info{atom_lines}[$idx];
        chomp($q_line);
        
        # Trim leading whitespace for cleaner output
        $q_line =~ s/^\s+//;

        printf "  [Index %3d] LAMMPS: ID=%-4d Type=%-2d | QE: %s\n", $idx, $l_id, $l_type, $q_line;
    }
    print "------------------------------------------------------\n";

    my %deleted_lookup = map { $_ => 1 } @to_delete;

    # --- Write New Files ---
    my $out_data = "$target_dir/${folder_name}.data";
    my $out_qe   = "$target_dir/${folder_name}.in";
    
    # 1. Write LAMMPS Data
    open(my $fh_d, '>', $out_data) or die "Cannot write $out_data: $!";
    
    my $new_total = $lammps_info{total_atoms} - $num;
    foreach my $line (@{$lammps_info{header}}) {
        if ($line =~ /(\d+)\s+atoms/) {
            print $fh_d "$new_total atoms\n";
        } else {
            print $fh_d $line;
        }
    }
    
    my $new_id = 1;
    for (my $idx = 0; $idx < $lammps_info{total_atoms}; $idx++) {
        next if $deleted_lookup{$idx}; 

        my $line = $lammps_info{atoms}[$idx]->[5];
        $line =~ s/^\s*(\d+)/$new_id/;
        print $fh_d $line;
        $new_id++;
    }
    close $fh_d;

    # 2. Write QE Input
    open(my $fh_q, '>', $out_qe) or die "Cannot write $out_qe: $!";
    
    foreach my $line (@{$qe_info{pre_atoms}}) {
        if ($line =~ /nat\s*=\s*\d+/) {
            $line =~ s/nat\s*=\s*\d+/nat = $new_total/;
        }
        print $fh_q $line;
    }
    
    for (my $idx = 0; $idx < $qe_info{total_atoms}; $idx++) {
        next if $deleted_lookup{$idx};
        print $fh_q $qe_info{atom_lines}[$idx];
    }
    
    print $fh_q $_ for @{$qe_info{post_atoms}};
    close $fh_q;

    print "Files created successfully.\n";
}

print "\nAll tasks completed.\n";

# --- Subroutines ---

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
            if ($line =~ /^\s*(\d+).+#\s*(\w+)/) { $elem_map{$2} = $1; }
        } elsif ($section eq "atoms_pre") {
            # Check for data start before pushing to header
            if ($line =~ /^\s*\d+/) { 
                $section = "atoms_body"; 
                redo; 
            }
            push @header, $line;
        } elsif ($section eq "atoms_body") {
            if ($line =~ /^\s*(\d+)\s+(\d+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)/) {
                # Store ID ($1) and Type ($2) explicitly for printing later
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
            if ($line =~ /K_POINTS|CELL_PARAMETERS|ATOMIC_SPECIES|^\s*[\&\/]/) {
                $section = "post";
                push @post, $line;
            } else {
                if ($line =~ /\w/) { push @atoms, $line; } 
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
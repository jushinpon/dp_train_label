#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(shuffle);
use File::Basename;
use Cwd 'abs_path';

# ==============================================================================
# USER CONFIGURATION
# ==============================================================================

# 1. Define arrays for your input files. 
#    Ensure the order matches: index 0 of data must match index 0 of QE input.

`rm -rf def-*`;#if needed, remove all old defect folders


my @data_sources = `find -L /home/jsp1/AlP/dp_train_new/initial -type f -name "*.data"|grep -E 'T50-P0|alp_2d-T0010'|grep -v _mp-`;#find all data files
map { s/^\s+|\s+$//g; } @data_sources;
@data_sources = sort @data_sources;
die "No data files to collect!\n" unless(@data_sources);

#my @data_sources = (
#    "/home/jsp1/AlP/dp_train_new/initial/Z_AlPNT_6x0-T50-P0/Z_AlPNT_6x0-T50-P0.data",
#    #"/home/jsp1/AnotherSystem/System2.data",
#    # Add more paths here...
#);

my @qe_sources = `find -L /home/jsp1/AlP/dp_train_new/initial -type f -name "*.in"|grep -E 'T50-P0|alp_2d-T0010'|grep -v _mp-`;#find all data files
map { s/^\s+|\s+$//g; } @qe_sources;
@qe_sources = sort @qe_sources;
die "No QE input files to collect!\n" unless(@qe_sources);

#my @qe_sources = (
#    "/home/jsp1/AlP/dp_train_new/initial/Z_AlPNT_6x0-T50-P0/Z_AlPNT_6x0-T50-P0.in",
#    #"/home/jsp1/AnotherSystem/System2.in",
#    # Add more paths here...
#);

# 2. Define defect parameters (applied to ALL files above)
my @delete_element = ("Al", "P"); # Element type to be deleted
my @delete_num     = (1, 1);              # Corresponding number to be deleted

# ==============================================================================
# MAIN PROCESSING LOOP
# ==============================================================================

# Verify array lengths match
if (scalar(@data_sources) != scalar(@qe_sources)) {
    die "Error: The number of data files (" . scalar(@data_sources) . ") does not match the number of QE input files (" . scalar(@qe_sources) . ").\n";
}

# Loop through each file pair
for (my $file_idx = 0; $file_idx < scalar(@data_sources); $file_idx++) {
    
    my $current_data = $data_sources[$file_idx];
    my $current_qe   = $qe_sources[$file_idx];

    # Resolve absolute paths
    my $abs_data_sour = abs_path($current_data);
    my $abs_qe_source = abs_path($current_qe);

    print "\n======================================================\n";
    print "PROCESSING FILE PAIR " . ($file_idx + 1) . "\n";
    print "Data: $current_data\n";
    print "QE:   $current_qe\n";
    print "======================================================\n";

    unless (-e $abs_data_sour && -e $abs_qe_source) {
        warn "Error: One or more input files not found for index $file_idx. Skipping.\n";
        next;
    }

    my ($base_name, $dir, $suffix) = fileparse($abs_data_sour, qr/\.[^.]*/);

    # 1. Read and Parse LAMMPS Data File
    my %lammps_info = read_lammps_data($abs_data_sour);

    # 2. Read and Parse QE Input File
    my %qe_info = read_qe_input($abs_qe_source);

    # Check consistency
    if ($lammps_info{total_atoms} != $qe_info{total_atoms}) {
        warn "Error: Atom count mismatch for $base_name! LAMMPS ($lammps_info{total_atoms}) vs QE ($qe_info{total_atoms}). Skipping.\n";
        next;
    }

    # 3. Apply Defects
    for (my $i = 0; $i < scalar(@delete_element); $i++) {
        my $elem = $delete_element[$i];
        my $num  = $delete_num[$i];
        
        my $prefix = "def-${elem}${num}";
        my $folder_name = "${prefix}-${base_name}";
        
        # Determine output directory (created in current working dir)
        my $current_work_dir = `pwd`;
        chomp($current_work_dir);
        my $target_dir = "$current_work_dir/$folder_name";
        
        print "\n  --- Task $i: Deleting $num atom(s) of element '$elem' ---\n";
        print "  Output: $target_dir\n";

        if (-d $target_dir) {
            `rm -rf "$target_dir"`; 
        }
        `mkdir -p "$target_dir"`;

        my $type_id = $lammps_info{element_map}{$elem};
        unless (defined $type_id) {
            warn "  Warning: Element '$elem' not found in LAMMPS Masses for $base_name. Skipping task.\n";
            next;
        }

        # Find candidates
        my @candidate_indices;
        for (my $idx = 0; $idx < $lammps_info{total_atoms}; $idx++) {
            if ($lammps_info{atoms}[$idx]->[1] == $type_id) {
                push @candidate_indices, $idx;
            }
        }

        if (scalar(@candidate_indices) < $num) {
            warn "  Warning: Not enough atoms to delete. Skipping task.\n";
            next;
        }

        my @shuffled = shuffle(@candidate_indices);
        my @to_delete = @shuffled[0 .. $num-1];
        
        # --- PRINT DELETED ATOM INFO ---
        print "  --------------------------------------------------\n";
        print "  Deleted Atoms:\n";
        foreach my $idx (@to_delete) {
            my $l_id   = $lammps_info{atoms}[$idx]->[0];
            my $l_type = $lammps_info{atoms}[$idx]->[1];
            my $q_line = $qe_info{atom_lines}[$idx];
            chomp($q_line);
            $q_line =~ s/^\s+//;
            printf "    [Idx %3d] LAMMPS ID=%-4d Type=%-2d | QE: %s\n", $idx, $l_id, $l_type, $q_line;
        }
        print "  --------------------------------------------------\n";

        my %deleted_lookup = map { $_ => 1 } @to_delete;

        # --- Write New Files ---
        my $out_data = "$target_dir/${folder_name}.data";
        my $out_qe   = "$target_dir/${folder_name}.in";
        
        # Write LAMMPS Data
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

        # Write QE Input
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

        print "  Created files successfully.\n";
    }
}

print "\nAll files processed.\n";

# ==============================================================================
# SUBROUTINES
# ==============================================================================

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
            if ($line =~ /^\s*\d+/) { 
                $section = "atoms_body"; 
                redo; 
            }
            push @header, $line;
        } elsif ($section eq "atoms_body") {
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
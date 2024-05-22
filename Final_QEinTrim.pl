#Kpoint sets 1 for z dimension
use strict;
use warnings;
use Cwd;
use POSIX;
#use lib '.';
#use elements;
#unlink "QEbyMatCld_summary.txt";
#open(my $DATA, ">QEbyMatCld_summary.txt");
#
#unlink "MatCld_crawling.txt";
#open(my $DATA1, ">MatCld_crawling.txt");

my $currentPath = getcwd();
my @QEin_MC = `find $currentPath/QEinByMatCld -name "*.in"`;
chomp @QEin_MC;

`rm -rf QE_trimmed4scf`;
`mkdir QE_trimmed4scf`;
#print $DATA1 "###Crawling QE input from Materials Cloud for:\n\n";
my $total = @QEin_MC;
my $counter = 0;
for my $qe (@QEin_MC){
    $counter++;
    my @magnet = `grep starting_magnetization $qe`;
    map { s/^\s+|\s+$//g; } @magnet;
    my $ntype = @magnet;
    my @pot = `grep -v '^[[:space:]]*\$' $qe|grep -A $ntype ATOMIC_SPECIES|grep -v ATOMIC_SPECIES|grep -v -- '--'`;
    map { s/^\s+|\s+$//g; } @pot;#need change $NF by the original QE input
    
    my $QEin_name = `basename $qe`;
    $QEin_name =~ s/^\s+|\s+$//g;
    my @pot_template = `grep -v '^[[:space:]]*\$' data2QE4MatCld/$QEin_name|grep -A $ntype ATOMIC_SPECIES|grep -v ATOMIC_SPECIES|grep -v -- '--'`;
    map { s/^\s+|\s+$//g; } @pot_template;#need change $NF by the original QE input
    my %ele_pot;
    for my $i (@pot_template){
        my @temp = split(/\s+/,$i);
        $ele_pot{$temp[0]} = "$temp[1] $temp[2]";
    }
   # for (keys %ele_pot){
   #     print "$_: $ele_pot{$_}\n";
   # }
   # die;
    my $kpoint_temp = `grep -v '^[[:space:]]*\$' $qe|grep -A 1 K_POINTS|grep -v K_POINTS`;
    $kpoint_temp =~ s/^\s+|\s+$//g;
    my @k_temp = split(/\s+/,$kpoint_temp);
    $k_temp[2] = 1;
    my $kpoint = join(" ",@k_temp);
    #print "$counter: $qe\n";
    #print "@magnet\n";
    #print "@pot\n";
    #print "$kpoint\n\n";
    #read the corresponding template    
    open my $in ,"< data2QE4MatCld/$QEin_name" or die "No data2QE4MatCld/$QEin_name";      
    my @QE_template =<$in>;
    close $in;
    map { s/^\s+|\s+$//g; } @QE_template;

    # find the key lines to modify
    my $pl;#id number with ATOMIC_SPECIES
    my $kl;#id number with K_POINTS
    my @start_mag;#ids with starting_magnetization

    for my $id (0..$#QE_template){
        if($QE_template[$id] =~ /ATOMIC_SPECIES/){
            $pl = $id;
        }
        elsif($QE_template[$id] =~ /K_POINTS/){
            $kl = $id;
        }
        elsif($QE_template[$id] =~ /starting_magnetization/){
            push @start_mag,$id;
        }
    }
    #begin trimming
    #pot
    for my $i (0 .. $ntype -1){
        my @temp = split(/\s+/,$pot[$i]);
        $QE_template[$pl+$i+1] = "$temp[0]     $ele_pot{$temp[0]}";
       # print "$pot[$i]\n";
       # print "$QE_template[$pl+$i+1]\n";
    }
    #kpoint
    $QE_template[$kl + 1] = $kpoint;
    #mag
     for my $i (0 .. $ntype -1){
        my $temp = $start_mag[$i];
        $QE_template[$temp] = $magnet[$i];
    }
    my $trimmed = join("\n",@QE_template);
    chomp $trimmed;
    open(FH, ">QE_trimmed4scf/$QEin_name" ) or die $!;
    print FH $trimmed;
    close(FH);
   # #$QEin_name =~ s/\.in//g;
   # chomp ($QEin_path, $QEin_name);
   # my $output = "QEinByMatCld/$QEin_name";#mv ...
   # unlink "$output";
   # unlink "input.in";
   # unlink "output.in";
   # `cp $qe input.in`; 
   # system("python ./QEinputByMatCld.py"); 
   # if($!){ print $DATA "$qe checked by Materials Cloud fail!\n";}
   # `mv output.in $output`;
}
#unlink "input.in";
#unlink "output.in";
#
#print "\n***Pleae check QEbyMatCld_summary.txt or the following content of QEbyMatCld_summary.txt:\n";
#print "\n***printing QEbyMatCld_summary.txt!!!\n";
#system("cat ./QEbyMatCld_summary.txt");
#print "\n\n*** If the above is empty, All QEinput files are good using Materials Cloud.\n";

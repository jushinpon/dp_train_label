##This module is developed by Prof. Shin-Pon JU at NSYSU on March 28 2021
package elements; 

use strict;
use warnings;

our (%element); # density (g/cm3), arrangement, mass, lat a , lat c

#$element{"N"} = [1.251,"hcp",14.007,3.861,6.265]; 
#$element{"Fe"} = [7.874,"bcc",55.845,2.8665,2.8665]; 
#$element{"Na"} = [0.968,"bcc",22.98977,4.2906,4.2906]; 
#$element{"O"} = [1.429,"fcc",15.9994,5.403,5.086]; 
#$element{"P"} = [1.823,"bcc",30.973761,1.25384,1.24896]; 


$element{"H"} = [0.0899, "hcp", 1.00794];
$element{"He"} = [0.1785, "hex", 4.002602];
$element{"Li"} = [0.535, "bcc", 6.941];
$element{"Be"} = [1.848, "hcp", 9.012182];
$element{"B"} = [2.34, "rhombohedral", 10.81];
$element{"C"} = [2.267, "hex", 12.011];
$element{"N"} = [0.00125, "hex", 14.007];
$element{"O"} = [0.00143, "cubic", 15.999];
$element{"F"} = [0.0017, "cubic", 18.998403163];
$element{"Ne"} = [0.9, "fcc", 20.1797];
$element{"Na"} = [0.971, "bcc", 22.98976928];
$element{"Mg"} = [1.738, "hcp", 24.305];
$element{"Al"} = [2.7, "fcc", 26.9815385];
$element{"Si"} = [2.33, "diamond", 28.085];
$element{"P"} = [1.82, "simple trigonal", 30.973761998];
$element{"S"} = [2.07, "orthorhombic", 32.06];
$element{"Cl"} = [0.003214, "orthorhombic", 35.45];
$element{"Ar"} = [0.001784, "fcc", 39.948];
$element{"K"} = [0.856, "bcc", 39.0983];
$element{"Ca"} = [1.54, "fcc", 40.078];
$element{"Sc"} = [2.99, "hcp", 44.955908];
$element{"Ti"} = [4.54, "hcp", 47.867];
$element{"V"} = [6.11, "bcc", 50.9415];
$element{"Cr"} = [7.15, "bcc", 52.];
$element{"Mn"} = [7.44, "bcc", 54.938044];
$element{"Fe"} = [7.87, "bcc", 55.845];
$element{"Co"} = [8.9, "hcp", 58.933194];
$element{"Ni"} = [8.908, "fcc", 58.6934];
$element{"Cu"} = [8.96, "fcc", 63.546];
$element{"Zn"} = [7.14, "hcp", 65.38];
$element{"Ga"} = [5.91, "orthorhombic", 69.723];
$element{"Ge"} = [5.323, "diamond", 72.630];
$element{"As"} = [5.72, "simple trigonal", 74.921595];
$element{"Se"} = [4.79, "hex", 78.96];
$element{"Br"} = [0.00759, "orthorhombic", 79.904];
$element{"Kr"} = [0.00375, "fcc", 83.798];
$element{"Rb"} = [1.532, "bcc", 85.4678];
$element{"Sr"} = [2.64, "fcc", 87.62];
$element{"Y"} = [4.47, "hcp", 88.90584];
$element{"Zr"} = [6.52, "hcp", 91.224];
$element{"Nb"} = [8.57, "bcc", 92.90637];
$element{"Mo"} = [10.22, "bcc", 95.94];
$element{"Tc"} = [11.5, "hcp", 97.90721];
$element{"Ru"} = [12.41, "hcp", 101.07];
$element{"Rh"} = [12.41, "fcc", 102.90550];
$element{"Pd"} = [12.02, "fcc", 106.42];
$element{"Ag"} = [10.49, "fcc", 107.8682];
$element{"Cd"} = [8.65, "hcp", 112.414];
$element{"In"} = [7.31, "tetragonal", 114.818];
$element{"Sn"} = [7.287, "diamond", 118.710];
$element{"Sb"} = [6.685, "simple trigonal", 121.760];
$element{"Te"} = [6.232, "hex", 127.60];
$element{"I"} = [4.93, "orthorhombic", 126.90447];
$element{"Xe"} = [0.005887, "fcc", 131.293];
$element{"Cs"} = [1.873, "bcc", 132.90545196];
$element{"Ba"} = [3.51, "bcc", 137.327];
$element{"La"} = [6.15, "hcp", 138.90547];
$element{"Ce"} = [6.77, "fcc", 140.116];
$element{"Pr"} = [6.77, "hcp", 140.90766];
$element{"Nd"} = [7.01, "hcp", 144.242];
$element{"Pm"} = [7.22, "hcp", 145.0];
$element{"Sm"} = [7.54, "hcp", 150.36];
$element{"Eu"} = [5.24, "bcc", 151.964];
$element{"Gd"} = [7.9, "hcp", 157.25];
$element{"Tb"} = [8.23, "hcp", 158.92535];
$element{"Dy"} = [8.55, "hcp", 162.500];
$element{"Ho"} = [8.8, "hcp", 164.93033];
$element{"Er"} = [9.07, "hcp", 167.259];
$element{"Tm"} = [9.32, "hcp", 168.93422];
$element{"Yb"} = [6.57, "fcc", 173.054];
$element{"Lu"} = [9.84, "hcp", 174.9668];
$element{"Hf"} = [13.31, "hcp", 178.49];
$element{"Ta"} = [16.65, "bcc", 180.94788];
$element{"W"} = [19.25, "bcc", 183.84];
$element{"Re"} = [21.02, "hcp", 186.207];
$element{"Os"} = [22.61, "hcp", 190.23];
$element{"Ir"} = [22.65, "fcc", 192.217];
$element{"Pt"} = [21.45, "fcc", 195.084];
$element{"Au"} = [19.32, "fcc", 196.966569];
$element{"Hg"} = [13.5336, "liquid", 200.592];
$element{"Tl"} = [11.85, "hcp", 204.38];
$element{"Pb"} = [11.34, "fcc", 207.2];
$element{"Bi"} = [9.747, "rhombohedral", 208.98040];
$element{"Th"} = [11.72, "hcp", 232.0377];
$element{"Pa"} = [15.37, "orthorhombic", 231.03588];
$element{"U"} = [18.95, "orthorhombic", 238.02891];
$element{"Np"} = [20.45, "orthorhombic", 237.0482];
$element{"Pu"} = [19.84, "monoclinic", 244.0642];
$element{"Am"} = [13.69, "double hexagonal close-packed", 243.0614];
$element{"Cm"} = [13.51, "double hexagonal close-packed", 247.0703];
$element{"Bk"} = [14.0, "double hexagonal close-packed", 247.0703];
$element{"Cf"} = [10.97, "double hexagonal close-packed", 251.0796];
$element{"Es"} = [8.84, "double hexagonal close-packed", 252.0829];

sub eleObj {# return properties of an element
    my $elem = shift @_;
    #print "\$elem: $elem\n";
    #print @{$element{"$elem"}}."\n";
    #print $element{"$elem"}->[0]."\n";
    #print $element{"$elem"}->[1]."\n";
    #print $element{"$elem"}->[2]."\n";
    return (\@{$element{"$elem"}});#return array reference only
}

1;               # Loaded successfully
#sub eleObj {# return properties of an element
#   my $elem = shift @_;
#   if(exists $element{"$elem"}){
#    return (@{$element{"$elem"}});      
#   }
#   else{
#      die "element information of \"$_\" is not listed in elements.pm.",
#      " You need to add Al according to the format of density (g/cm3), arrangement, mass, lat a , lat c. ",
#      ' For example, $element{"Nb"} = [8.57,"bcc",92.90638,3.30,3.30]'."\n"; 
#   }
#}
#1;               # Loaded successfully

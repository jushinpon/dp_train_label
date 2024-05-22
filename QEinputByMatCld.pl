use strict;
use warnings;
use Cwd;
use POSIX;
#use lib '.';
#use elements;
unlink "QEbyMatCld_summary.txt";
open(my $DATA, ">QEbyMatCld_summary.txt");

unlink "MatCld_crawling.txt";
open(my $DATA1, ">MatCld_crawling.txt");

my $currentPath = getcwd();
my @QEin = `find $currentPath/data2QE4MatCld -name "*.in"`;
chomp @QEin;

`rm -rf QEinByMatCld`;
`mkdir QEinByMatCld`;
print $DATA1 "###Crawling QE input from Materials Cloud for:\n\n";
my $total = @QEin;
my $counter = 0;
for my $qe (@QEin){
    $counter++;
    print $DATA1 "$counter: $qe -> $counter/$total\n";
    print "Crawling $qe status: $counter/$total\n";
    my $QEin_path = `dirname $qe`;
    my $QEin_name = `basename $qe`;
    #$QEin_name =~ s/\.in//g;
    chomp ($QEin_path, $QEin_name);
    my $output = "QEinByMatCld/$QEin_name";#mv ...
    unlink "$output";
    unlink "input.in";
    unlink "output.in";
    `cp $qe input.in`; 
    system("python ./QEinputByMatCld.py"); 
    if($!){ print $DATA "$qe checked by Materials Cloud fail!\n";}
    `mv output.in QEinByMatCld/$QEin_name`;
}
unlink "input.in";
unlink "output.in";

print "\n***Pleae check QEbyMatCld_summary.txt or the following content of QEbyMatCld_summary.txt:\n";
print "\n***printing QEbyMatCld_summary.txt!!!\n";
system("cat ./QEbyMatCld_summary.txt");
print "\n\n*** If the above is empty, All QEinput files are good using Materials Cloud.\n";

use strict;
use warnings;
open (IN, 'fs21Ucounts.csv') or die $!;
my $skip = 0;
while (<IN>){
    #      Name             Sequence      Embryo    mL1     mL2       mL3       mL4       yAdult    yAdult
    # † 21UR-1    TGGTACGTACGTTAACCGTGC         0       0        0         0          0        0         0            0
    #   21UR-2    TGGGAAATTCGAATAATATAT         0       0        0         1          2        1         1            5
    $skip=1 if ($_ =~/Name/);
    next unless $skip;
    chomp $_;
    my @arr = split /\s+/,$_;
    for (@arr) { shift if s/\s+|\*|†//g }
    @arr = grep { $_ !~ /^\s?$/ } @arr;
    print join "\t",@arr;
    print "\n";
}

close IN;

use strict;
use warnings;
open (IN, '/home/philip/other_data/gff/frank_slack_may_2009_data/fs21Ucounts.clean.csv') or die $!;
my %fs21Us;
while (<IN>){
    # Name	Sequence	Embryo	mL1	mL2	mL3	mL4	yAdult	yAdult
    # 21UR-1	TGGTACGTACGTTAACCGTGC	0	0	0	0	0	0	0	0
    next if ($_ =~/Name/);
    chomp $_;
    my @ar = split /\t/,$_;
    $fs21Us{uc($ar[0])}=1;
}
close IN;

print scalar keys %fs21Us,"\n";

open (IN, '/home/philip/other_data/gff/WB_21U.gff') or die $!;
while (<IN>){
    $_ =~ m/.*Locus \"(.*)\"/;
    my $u = uc($1);
    if (!$fs21Us{$u}){
	print $u,"\n";
    }
}
close IN;

use strict;
use warnings;

# We don't know what 21Us are correlated with what Tx's
# 
# We want to be able to correlate total number of 21U molecules to
# some measure of Expression data
# 
# method : 
# look at all expression data, take avg. 
# pick some number; rad, that is the max distance from avg
# count all guys that lie outside that distance 
# for 0 - XX counted 21Us
# plot num 21Us v num guys outside avg with increasing 21U counts.


# look at all expression data, take avg. 
my $geneToLogExprn = getLogExprn();

my $avgLogExprsn = getAvgExprn($geneToLogExprn);

my ($U21CountsYoungAdult, $U21CountsTotal) = get21UCounts();

my $U21FromGFF = load21InfoFromGFF();

foreach (keys %{$U21FromGFF}){
    my %U21= %{$U21FromGFF->{$_}};
    foreach (keys %U21){
	print $U21{
    
#    push (@a, $_) if $U21CountsYoungAdult->{$_}==0
}


sub load21InfoFromGFF{

    open IN, '</mnt/disk1/philip_scratch/genomes/celegans/annot/WB_21U.gff' or die $!;
    # CHROMOSOME_I	gene	gene	31523	31543	.	+	.	Gene "WBGene00170953" ; Interpolated_map_position "-21.244" ; Locus "21ur-15479"
    my %h;
    my ($prevStart, $prevStrand, $prevUname);
    while (<IN>){
	chomp $_;
	my @a = split "\t",$_;
	next unless $a[1] eq 'gene';
	my $uName = $a[$#a];
	if ($uName =~ m/Locus "(.*)"/){
	    $uName = $1;
	}
	$a[0] =~ s/CHROMOSOME_//;
	
	my $strand = $a[6];
	$h{$uName}{csome}=$a[0];
	$h{$uName}{strand} = $strand;

	my $start;
	$strand eq '+' ? $start = $a[3]: $start = $a[4];
	$h{$uName}{start} = $start;

	# first time we run loop or csome or strand changes
	if ( !(defined $prevStrand) || ($strand ne $prevStrand)){
	    undef $prevStart;
	}

	if (defined $prevStart)
	{
	    $h{$uName}{nearest}=$prevUname;
	    my $thisGap = abs($prevStart-$start);
	    $h{$uName}{gap}=$thisGap;

	    my $prevGap = $h{$prevUname}{gap};
	    if ($thisGap < $prevGap){
		$h{$prevUname}{nearest}=$uName;
		$h{$prevUname}{gap}=$prevGap;
	    }
	    $prevGap = $thisGap;
	}

	if (!defined $h{$uName}{nearest}){
	    $h{$uName}{nearest}=$prevUname;
	}

	$prevStrand = $strand;
	$prevStart = $start;
	$prevUname = $uName;
    }
    \%h;
}

sub countPopOutsideRad{
    my ($thisRad,$h)=@_;
}    

sub get21UCounts{
    my (%h, %i);
    open (IN, '</mnt/disk1/philip_scratch/other_data/frank_slack_may_2009_data/fs21Ucounts.clean.csv') or die $!;
    while (<IN>){
	chomp$_;
	next if $_ =~ /^Name|^Total/;
	my @a = split /\t/,$_;
	$h{$a[0]}=$a[6]+$a[7]+$a[8];
	$i{$a[0]}=$a[9];
    }
    close IN;
    return (\%h,\%i);
}

sub getAvgExprn{
    my $a;
    my $toAvg=$_[0];
    foreach (keys %{$toAvg}){
	$a += $toAvg->{$_}
    }
    return sprintf("%.3f", $a/(scalar keys %{$toAvg}));
}

sub getLogExprn{
    my %geneToLogExprnL;
    open (MA, '</home/philip/other_data/prg1_microarray.clean.50.txt') or die $!;
    while (<MA>){	
	chomp $_;
	$_ =~ s/"//g;
	my @a = split /\s/,$_;
	my $geneName = $a[1];
	my $lfc = $a[7];
	if ($geneToLogExprnL{$geneName}){
	    $geneToLogExprnL{$geneName} = ($geneToLogExprnL{$geneName}+$lfc)/2;
	}
	else {
	    $geneToLogExprnL{$geneName}=$lfc;
	}
    }
    close MA;
    return \%geneToLogExprnL;
}

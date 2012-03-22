use strict;
use warnings;

# loads of sucking hacky shit in here... 
# load up all Reinke stuff. 
# map geneName to log(exp rat)

my %geneToExp;
my %cdsToExp; 
my $lowVal = -99999;
open (IN, "/home/philip/other_data/prg1_alldata.clean.txt") or die $!;
while  (<IN>)
{
    # "Name" "Gene name without alt splice info" "Gene" "Chromosome" "Strand" "Start" "Stop" "average ratio"
    # "D2030.6" "D2030.6" "prg-1" "I" -1 7585320 7588362 3.91208408408408
    #  ^^CSD 0   ^^gene 1  2       3   4  5      6       7
    chomp $_;
    $_ =~ s/"//g;
    next if $_ =~/^Name/;
    my (undef,$geneName,undef,undef,undef,undef,undef,$ratio) = split /\s/,$_;
    my $logRatio;
    if ($geneToExp{$geneName})	# if we have Reinke data for this gene
    {
	my $l = log_2($ratio);	# look up log of the ratio
	if ($geneToExp{$geneName} == $lowVal){
	    $logRatio = log_2($ratio)
	}
	elsif ($l != $lowVal) {
	    $logRatio = ($l+$geneToExp{$geneName})/2;
	}
	else {
	    $logRatio = $geneToExp{$geneName};
	}
    }
    else {
	$logRatio = log_2($ratio)
    }
  
    if ($logRatio == $lowVal) {
	$geneToExp{$geneName} = $logRatio;
    }
    else {
	$geneToExp{$geneName} = sprintf("%.4f", $logRatio);
    }
}
close IN or warn $!;

my @bowtieFiles = `ls /home/philip/other_data/bowtie/bin_offrate_1/outputs/*out`;

foreach my $fileName (@bowtieFiles)
{
    chomp $fileName;
    open IN, $fileName or die $!;

    # /home/philip/other_data/Cel_21U_v_3_prime0MM_nofw.out
    my @fileNames = split /\//,$fileName;
    my $outFileName = $fileNames[$#fileNames].'_Reinke_expression';
    my $noExpOutFileName = $fileNames[$#fileNames].'_no_exprsn_avail';
    open OUT, ">/home/philip/other_data/bowtie/bin_offrate_1/cross_refd_outputs/$outFileName" or die $!;
    open NO_EXP_OUT, ">/home/philip/other_data/bowtie/bin_offrate_1/cross_refd_outputs/$noExpOutFileName" or die $!;


    my %URNAs;
    while (<IN>)
    {
        # 5 prime
        # 21ur-13748|F44F4.14|II|10886034|10886054|+	-	WBGene00003687|H27C11.1a|H27C11.1|H27C11.1a	130

	# CDS
        # 21ur-14112|Y71G12B.34|I|1653767|1653787|+	-	F55D1.1	196

	# 3 prime
	# 21ur-14201|F21C3.9|I|7286319|7286299|-	-	ContigIV.elegans.jigsaw.1595.mRNA1595_g11275.t2|C28C12.11a_C28C12.11b|+	745

	chomp $_;
	my ($U21, $strand, $alignedTo, $offset)= split /\t/,$_;

	my @alignedToArr = split /\|/, $alignedTo;

	my $geneName;
	if ($alignedTo !~ /\|/) {
	    $geneName = $alignedTo
	}
	elsif ($alignedTo =~ /^WB/) { 
	    $geneName = $alignedToArr[2]
	}
	else {
	    $geneName = toGeneName($alignedToArr[1])
	}

	if ($geneName eq'') {
	    print "Cannot look up gene name for $alignedTo \n";
	    exit;
	}

	my %h;

	$h{geneHit}= $geneName;

	if (!defined $geneToExp{$geneName}){
	    print NO_EXP_OUT join "\t",($_,"\n");
	    next;
	}

	$h{logRatio} = $geneToExp{$geneName};

	my @U21arr = split /\|/,$U21;
	$h{name} = $U21arr[0];
	$URNAs{$geneName} = \%h;
    }
    close IN or warn $!;

    foreach my $l (sort {$URNAs{$a}->{logRatio} <=> $URNAs{$b}->{logRatio}} keys %URNAs)
    {
	my $href = $URNAs{$l};
	if (!$href->{geneHit})
	{
	    print NO_EXP_OUT "Cannot look up a geneName for -- $l --\n";
	    next;
	}
	# my $logRat;
	# if (!$href->{logRatio}) { $logRat = 0 }
	# else {$logRat = }
	print OUT join "\t", ($href->{geneHit}, $href->{name}, $href->{logRatio},"\n");
    }
    close OUT or warn;
    close NO_EXP_OUT or warn;
}

sub log_2
{
    my $n = shift;
    return $lowVal if $n =~ /DIV/ || $n <= 0;
    return log($n)/log(2);
}

sub toGeneName{
    my $name = shift;
    if ($name =~m/-1/) { 
	return 0;
	# print "Exiting !! ---> $_ <---\n"; 
	# exit
    }

    return "" if !defined $name || $name eq '';
    my @a = split /_/,$name;
    return $name if (length($a[0]) == 1);
    if (!$a[0]) {print "nothing for $name!!! Exiting"; exit}
    my @b = split /\./,$a[0];
    if (!$b[1]) {print "nothing for $name!!! Exiting"; exit}
    $name = join '.',@b[0..1];
    $name =~ s/[a-z]$//;
    return $name;
}
    
######

use strict;
use warnings;

# loads of sucking hacky shit in here... 
# load up all Reinke stuff. 
# map geneName to log(exp rat)

my %geneToExp;
my %cdsToExp; 

open (IN, "/home/philip/other_data/prg1_alldata.clean.txt") or die $!;
while  (<IN>)
{
    # "Name" "Gene name without alt splice info" "Gene" "Chromosome" "Strand" "Start" "Stop" "average ratio"
    # "D2030.6" "D2030.6" "prg-1" "I" -1 7585320 7588362 3.91208408408408
    #  ^^CSD 0   ^^gene 1  2       3   4  5      6       7
    chomp $_;
    $_ =~ s/"//g;
    next if $_ =~/^Name/;
    my @line = split /\s/,$_;
    my $logRatio = 0;
    if ($geneToExp{$line[1]})
    {
	my $l = log_2($line[7]);
	$logRatio = ($l+$geneToExp{$line[1]})/2;
    }
    else { $logRatio = log_2($line[7]) }
  
    
    $geneToExp{toGeneName($line[1])} = $logRatio;

}
close IN or warn $!;

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

my @bowtieFiles = `ls /home/philip/other_data/*out`;

foreach my $fileName (@bowtieFiles)
{
    chomp $fileName;
    open IN, $fileName or die $!;

    # /home/philip/other_data/Cel_21U_v_3_prime0MM_nofw.out
    my @fileNames = split /\//,$fileName;
    my $outFileName = $fileNames[$#fileNames].'_Reinke_expression';
    my $noExpOutFileName = $fileNames[$#fileNames].'_no_exprsn_avail';
    open OUT, ">/home/philip/other_data/$outFileName" or die $!;
    open NO_EXP_OUT, ">/home/philip/other_data/$noExpOutFileName" or die $!;


    my %URNAs;
    while (<IN>)
    {
	# 21ur-15425|F18C5.11|II|6570425|6570405|-	-	ContigII.elegans.jigsaw.1733.mRNA1733|F18C5.1|+	2	TTCCTCAATCCAGAAACATTA	IIIIIIIIIIIIIIIIIIIII	0
	#                                                                                             ^^ gene hit
	# 21ur-12947|F02C12.6|X|13397458|13397438|-	-	F02C12.3	37	AGCAGAAGCTCAAGTGGAGAA	IIIIIIIIIIIIIIIIIIIII	0
	# 21ur-11631|D1007.21|I|4584789|4584769|-	-	D1007.6.2|-1|4584516|4584960	253	GGTACAACGAAATTCTTCCGA	IIIIIIIIIIIIIIIIIIIII
	chomp $_;

	my @line = split /\|/,$_;

	my $geneName;

	if (scalar @line == 6) {
	    my @a = split /\t/,$line[5];
	    $geneName = toGeneName($a[2]);
	}
	elsif (scalar @line ==8 || scalar @line ==9) {
	    $geneName = toGeneName($line[6]);
	    if (!$geneName) {
		my @a = split /\t/,$line[5];
		$geneName = toGeneName($a[2]);
	    }
	}
	else {
	    print "exiting on $_\n\n"; exit
	}
	my %h;


	$h{geneHit}= $geneName;

	if (!defined $geneToExp{$geneName}){
	    print NO_EXP_OUT join "\t",(@line,"\n");
	    $h{logRatio}=0;
	}
	else { 
	    $h{logRatio} = $geneToExp{$geneName}
	}

	$h{name} = $line[0];

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
	my $logRat;
	if (!$href->{logRatio}) { $logRat = 0 }
	else {$logRat = $href->{logRatio}}
	print OUT join "\t", ($href->{geneHit}, $href->{name}, $logRat,"\n");
    }
    close OUT or warn;
    close NO_EXP_OUT or warn;
}

sub log_2
{
    my $n = shift;
    return -999999 if $n =~ /DIV/ || $n <= 0;
    return log($n)/log(2);
}

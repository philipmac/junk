use strict;
use warnings;

# all this does is take all of the outputs from Bowtie and make them more simple to play with on the CLI
# nothing more!

my @bowtieFiles = `ls /home/philip/other_data/bowtie/bin_offrate_1/outputs/*out`;
foreach my $fileName (@bowtieFiles)
{
    chomp $fileName;
    open IN, $fileName or die $!;

    # /home/philip/other_data/Cel_21U_v_3_prime0MM_nofw.out
    my @fileNames = split /\//,$fileName;
    my $outFileName = $fileNames[$#fileNames].'_parsed';

    open OUT, ">/home/philip/other_data/bowtie/bin_offrate_1/parsed_outputs/$outFileName" or die $!;

    my %URNAs;
    while (<IN>)
    {
        # 5 prime
        # 21ur-13748|F44F4.14|II|10886034|10886054|+	-	WBGene00003687|H27C11.1a|H27C11.1|H27C11.1a	130

	# CDS
        # 21ur-14112|Y71G12B.34|I|1653767|1653787|+	-	F55D1.1	196

	# 3 prime
	# 21ur-14201|F21C3.9|I|7286319|7286299|-	-	ContigIV.elegans.jigsaw.1595.mRNA1595_g11275.t2|C28C12.11a_C28C12.11b|+	745
	# 21ur-5849|VY10G11R.50|IV|16460711|16460691|-	-	WBsf218491|C35E7.5a_C35E7.5b|+	199

	chomp $_;
	my ($U21, $strand, $alignedTo, $offset)= split /\t/,$_;

	my @alignedToArr = split /\|/, $alignedTo;
	my @U21Arr = split /\|/, $U21;

	my $geneName;
	if ($alignedTo !~ /\|/) {
	    $geneName = $alignedTo;
	}
	elsif ($alignedTo =~ /^WBGene/) { 
	    $geneName = $alignedToArr[2];
	}
	else {
	    $geneName = toGeneName($alignedToArr[1]);
	}

	if ($geneName eq'') {
	    print "Cannot look up gene name for $alignedTo \n";
	    exit;
	}
	print OUT join "\t",($geneName,$U21Arr[0]);
	print OUT "\n";
    }
    close OUT;
}


sub toGeneName{
    my $name = shift;
    if ($name =~m/-1/) { 
#	return 0;
	print "Exiting !! ---> $_ <---\n"; 
	exit
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

use strict;
use Bio::DB::GFF;

my $cds = Bio::DB::GFF::Aggregator->new(-method=> 'pcds', 
					-sub_parts=>['three_prime_UTR:Coding_transcript',
						      'coding_exon:Coding_transcript',
						      'five_prime_UTR:Coding_transcript'
						     ]); # CDS:curated',
my $db      = Bio::DB::GFF->new(-adaptor => 'dbi::mysqlopt',
                                -dsn     => 'dbi:mysql:WS229', 
                                -aggregator=> $cds,
                                -user => 'root',
                                -password => 'argonaute'

);
my $longestPromotor = 999;

my $file = IO::File->new("> promotors.fa") or die $!;
foreach my $csome (qw /CHROMOSOME_I CHROMOSOME_II CHROMOSOME_III CHROMOSOME_IV CHROMOSOME_V CHROMOSOME_X/)
{
    my %occupiedPositions;
    my $segment = $db->segment($csome); # , 1=>100000
    $segment->absolute(1);
    
    my @transcripts = $segment->contained_features('pcds');
    foreach my $transcript (@transcripts)
    {
	foreach my $loc($transcript->abs_low..$transcript->abs_high)
        {
	    $occupiedPositions{$loc}=1;
	}
    }

    my %alreadyUsed;
    foreach my $transcript (@transcripts)
    {
	my ($thisPromoterStart,$thisPromoterEnd,$thisPromLen) = (0,0,0);
	$transcript->absolute(1);

	if ($transcript->strand == 1) { 
	    $thisPromoterStart=$transcript->abs_low;
	    $thisPromoterEnd = $thisPromoterStart-1;
	}
	else { 
	    $thisPromoterStart=$transcript->abs_high;
	    $thisPromoterEnd = $thisPromoterStart+1;
	}

	while(1)
	{
	    if (!$occupiedPositions{$thisPromoterEnd} && ($thisPromLen < $longestPromotor) && !$alreadyUsed{$thisPromoterEnd})
	    {
		$alreadyUsed{$thisPromoterEnd}=1;
		if ($transcript->strand ==1) { $thisPromoterEnd-- }
		else { $thisPromoterEnd++ }
		$thisPromLen++;
	    }
	    else {last}
	}
#	print join ",", ($transcript, $thisPromoterStart, $thisPromoterEnd,$thisPromLen,$transcript->strand,"\n");
	if ($transcript->strand ==1) { $thisPromoterEnd++ }
	else { $thisPromoterEnd-- }
#	print join ",", ($transcript,$thisPromoterStart,$thisPromoterEnd,$transcript->strand,"\n");
	next if $thisPromLen==0;
	my $promSeg = $segment->subseq($thisPromoterStart, $thisPromoterEnd);
	print $file join "|",(">".$transcript->name, $transcript->strand, $thisPromoterStart, $thisPromoterEnd);
	print $file "\n",$promSeg->dna,"\n";
    }
}
$file->close;

use strict;
use Bio::DB::GFF;
use Data::Dumper;

# UTRs seem not to be associated with any CDS explicitly.  (lame)  Still going to agg them because they might have introns.
my $cds_utr = Bio::DB::GFF::Aggregator->new(-method=> 'utr', 
					    -sub_parts=>[('three_prime_UTR')]); # CDS:curated',

my $db      = Bio::DB::GFF->new(-adaptor => 'dbi::mysqlopt',
				-dsn     => 'dbi:mysql:WS229', 
				-aggregator=> [$cds_utr],
				-user => 'root',
				-password => 'argonaute'

);

my $file = IO::File->new("> UTRomeUTRs.fa") or die $!;
my $junk_file = IO::File->new("> bad_utr.txt") or die $!;

print "Enter\n";

foreach my $csome (qw /CHROMOSOME_I CHROMOSOME_II CHROMOSOME_III CHROMOSOME_IV CHROMOSOME_V CHROMOSOME_X/) # 
{
    my (%spaceUsedAlreadyPlus, %spaceUsedAlreadyMinus);
    
    my $segment = $db->segment($csome); #
    $segment->absolute(1);
    my (%posToBasePos, %posToBaseNeg);

    # zero based or 1 based? <<-- seemingly 1
    print "getting dna\n";
    my $dnaStr = $segment->dna;
    my @dna = split //,$dnaStr;
    $dnaStr = '';

    my $revcom =  join('',@dna);  #reverse join('',@dna); reverse join('',@dna);
    $revcom =~ tr/ACGTacgt/TGCAtgca/;
    my @revDNA = split //,$revcom;
#    print @revDNA[0..100];
#    exit;
    print "Done DNA\n";

    print "getting utr\n";
    my @utrs = $segment->contained_features('utr'); # arg agg above
    print "Done UTR, looking up info\n";

    foreach my $utr (@utrs)
    {
	# print "-->", $utr->sub_SeqFeature;
	# get the strand, ignore stuff that's on the other strand.
	$utr->absolute(1);	# abs_start being used below
	my $utr_strand = $utr->strand;
	
	my $assocCDSs_href = lookup_CDS($utr, $segment);

	foreach my $part ($utr->sub_SeqFeature)
        {	    
#	    print $part->abs_low,"..",$part->abs_high,"\n";
	    foreach my $loc ($part->abs_low..$part->abs_high)
	    {
		my %h;
		$h{cds} = $assocCDSs_href;
		$h{utr} = {$utr->name=>1};
		$h{strand} = $utr_strand;

		#print $loc,"->",$dna[$loc-1],"\n";
		if($utr_strand == 1)
		{
		    $h{nuc} = $dna[$loc-1];
		    $posToBasePos{$loc} = \%h;
		}
		else
		{
		    $h{nuc}=$revDNA[$loc];
		    $posToBaseNeg{$loc} = \%h;
		}
	    }
	}
	# print $junk_file $utr->name." Space used twice + \n" if ($spaceUsedAlreadyPlus{$_});
	# print $junk_file $utr->name." Space used twice - \n" if ($spaceUsedAlreadyMinus{$_});
    }
    
    print "done info\n";
    @revDNA = @dna = @utrs = ();


    my $catDNA;
    my (%cdsToUniq, %utrToUniq);

#, @posToBaseNeg
    my $prevLoc = 0;
    foreach my $loc (sort {$a<=>$b} keys %posToBasePos)
    {

	#print $loc,"\n";
#c489819,t489820,c489821,g489822,
	if ($prevLoc+1 != $loc) 	# if we hit an elt in the arr that has no object attached, jumped off the end of a UTR, need to print out and dump. or last itration
        {
	    my $cdsStr = join "_",keys %cdsToUniq;
	    my $utrStr = join "_",keys %utrToUniq;
	    if ($cdsStr ne '')
	    {
		print $file join '|', (">$utrStr",$cdsStr,'+');
		print $file "\n",$catDNA,"\n";
	    }

	    %cdsToUniq = ();
	    %utrToUniq = ();
	    $catDNA ='';
	}

	$prevLoc=$loc;

	my %obj=%{$posToBasePos{$loc}};
	# my $d = Data::Dumper->new([$obj]);
	# print $d->Dump;
	# next;
	
	foreach(keys %{$obj{cds}})
	{
	    $cdsToUniq{$_}=1;
	}

	foreach(keys %{$obj{utr}})
	{
	    $utrToUniq{$_}=1;
	}
	$catDNA .= $obj{'nuc'};
#	print $loc,",",$obj{'nuc'},"\n";
	    #('>'.$utr->name, $csome, $utr->abs_start, $utr->abs_end, $strand);
	    #@cds=keys(
    }
    
    my $cdsStrP = join "_",keys %cdsToUniq;
    my $utrStrP = join "_",keys %utrToUniq;
    print $file join '|', (">$utrStrP",$cdsStrP,'+');
    print $file "\n",$catDNA,"\n";

    %cdsToUniq = ();
    %utrToUniq = ();
    $catDNA ='';

    foreach my $loc (sort {$b<=>$a} keys %posToBaseNeg)
    {
	if ($prevLoc-1 != $loc)
        {
	    my $cdsStr = join "_",keys %cdsToUniq;
	    my $utrStr = join "_",keys %utrToUniq;
	    if ($cdsStr ne '')
	    {
		print $file join '|', (">$utrStr",$cdsStr,'-');
		print $file "\n",$catDNA,"\n";
	    }
	    %cdsToUniq = ();
	    %utrToUniq = ();
	    $catDNA ='';
	}

	$prevLoc=$loc;

	my %obj=%{$posToBaseNeg{$loc}};
	# my $d = Data::Dumper->new([$obj]);
	# print $d->Dump;
	# next;
	
	foreach(keys %{$obj{cds}})
	{
	    $cdsToUniq{$_}=1;
	}

	foreach(keys %{$obj{utr}})
	{
	    $utrToUniq{$_}=1;
	}
	$catDNA .= $obj{'nuc'};
#	print $loc,",",$obj{'nuc'};
    }
}


$file->close;			
$junk_file->close;

sub lookup_CDS() #look up the gene... # want to find a cds that has a stop near or at UTR start
{
    my ($utr, $segment) = @_;
    my $searchSeg = $segment->subseq($utr->abs_start-10, $utr->abs_end+10);
    $searchSeg->absolute(1);

    my @allCDS = $searchSeg->features('coding_exon:curated');

    my (%good_cds, %utr_to_cds);

    foreach my $cds (@allCDS)
    {
	next if $cds->strand != $utr->strand;
	if ($utr->strand == 1)
	{
	    if (($cds->abs_end >= $utr->abs_start-10) && ($utr->abs_end > $cds->abs_end))
	    {
		$utr_to_cds{$cds->name}=1;
	    }
	    
	}
	else
	{
	    my ($cds_end, undef ) = sort {$a<=>$b} ($cds->abs_start, $cds->abs_end);
	    my ($utr_stop, $utr_start ) = sort {$a<=>$b} ($utr->abs_start, $utr->abs_end);
#	    $cds_of_interest{$cds->name} = 1 if (($cds_end <= $utr_start+10) && ($utr_stop < $cds_end));  
	    $utr_to_cds{$cds->name}=1;
	}
    }

    if (!%utr_to_cds)	
    {  
	foreach my $part (@allCDS)
	{
	    print $junk_file $part->strand," ",$part->abs_start," ",$part->abs_end," ",$part->type," ",$part->name, "\n";
	}
	print $junk_file "utr seg->";
	print $junk_file join "..", ('I',$utr->abs_start, $utr->abs_stop, $utr->strand, "\n\n");
    }

    return \%utr_to_cds;
}	

		# if (! %{$posToBasePos[$loc]{'cds'}})
                # {
		#     %{$posToBasePos[$loc]{cds}) = {'cds'=>{""}};
		# }   
	    # my (undef, $cds_end) = sort {$a<=>$b} ($cds->abs_start, $cds->abs_end);
	    # my ($utr_start, $utr_stop) = sort {$a<=>$b} ($utr->abs_start, $utr->abs_end);
	    # $cds_names{$cds->name} = 1 if (($cds_end >= $utr_start-10) && ($utr_stop > $cds_end));

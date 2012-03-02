use strict;
use warnings;
use Bio::DB::GFF;
use Data::Dumper;

# UTRs seem not to be associated with any CDS explicitly.  (lame)  Still going to agg them because they might have introns.
my $cds_utr = Bio::DB::GFF::Aggregator->new(-method=> 'cds', 
                                            -sub_parts=>[('coding_exon:Coding_transcript')]);

my $db      = Bio::DB::GFF->new(-adaptor => 'dbi::mysqlopt',
                                -dsn     => 'dbi:mysql:WS229', 
                                -aggregator=> [$cds_utr],
                                -user => 'root',
                                -password => 'argonaute'

);

my $file = IO::File->new("> CDSSuperSet.fa") or die $!;

sub rev_com
{
    my $in = shift;

    my @dna = split //,$in;
    my $revcom = reverse join('',@dna);
    $revcom =~ tr/ACGTacgt/TGCAtgca/;
    return $revcom;
}

sub superset
{
    my %a = %{$_[0]};
    my %b = %{$_[1]};
    return my %h if (!keys %a && !keys %b);
    return \%b if !keys %a;
    return \%a if !keys %b;
    foreach (keys %a) 
    { 
	$b{$_} = 1
    }
    return \%b;
}

sub hashes_intersect
{
    my %a = %{$_[0]};
    my %b = %{$_[1]};
    if ((!keys %a) && (!keys %b))
    {	print "exiting!", keys %a, " ",keys %b;
	exit;
    }

    foreach (keys %a)
    {
	return 1 if $b{$_}
    }

    return 0;
}

sub printOut
{
    my %g = %{$_[0]};
    my $d = $_[1];
    my $topLine = join ('|', (keys %g));
    print $file '>'.$topLine."\n".$d."\n";
}

foreach my $csome (qw /CHROMOSOME_I CHROMOSOME_II CHROMOSOME_III CHROMOSOME_IV CHROMOSOME_V CHROMOSOME_X/) # 
{
    my @locMapperPlus;
    my @locMapperMinus;

    my $segment = $db->segment($csome);
    $segment->absolute(1);

    my $dnaStr = $segment->dna;
    my @dna = split //,$dnaStr;
    $dnaStr = '';

    $segment->absolute(1);
    my @allCDS = $segment->features('coding_exon:curated');
    foreach my $exon (@allCDS)
    {
	$exon->absolute(1);
	my $exonName = stripName($exon->name);

	foreach my $loc ($exon->abs_start..$exon->abs_end)
	{
	    if ($exon->strand eq '1') { $locMapperPlus[$loc]->{$exonName}=1 }
	    else { $locMapperMinus[$loc]->{$exonName}=1 } 
	}
    }

    my $dnaSt='';
    my $prevDNA='';
    my %prevGenes = ();
    my %genes = ();

    for my $i (0 .. $#locMapperPlus)
    {
	if (defined $locMapperPlus[$i])
	{
	    $dnaSt .= $dna[$i];
	    foreach (keys %{$locMapperPlus[$i]})
	    { 
		$genes{$_}=1;
	    }
	}
	elsif($dnaSt)
	{
	    if(hashes_intersect(\%prevGenes, \%genes) || (!keys %prevGenes))
	    {
		%prevGenes = %{superset(\%prevGenes, \%genes)};
		$prevDNA .= $dnaSt;
	    }
	    else
	    {
		printOut(\%prevGenes, $prevDNA);
		$prevDNA = $dnaSt;
		%prevGenes = %genes;
	    }
	    $dnaSt = '';
	    %genes = ();
	}
    }

    printOut(\%prevGenes, $prevDNA);
    $prevDNA = $dnaSt;
    %prevGenes = %genes;
    $dnaSt = '';
    %genes = ();

    for my $i (0 .. $#locMapperMinus)
    {
	if (defined $locMapperMinus[$i])
	{
	    $dnaSt .= $dna[$i];
	    foreach (keys %{$locMapperMinus[$i]})
	    { 
		$genes{$_}=1;
	    }
	}
	elsif($dnaSt)
	{
	    if(hashes_intersect(\%prevGenes, \%genes) || (!keys %prevGenes))
	    {
		%prevGenes = %{superset(\%prevGenes, \%genes)};
		$prevDNA .= $dnaSt;
	    }
	    else
	    {
		printOut(\%prevGenes, rev_com($prevDNA));
		$prevDNA = $dnaSt;
		%prevGenes = %genes;
	    }
	    $dnaSt = '';
	    %genes = ();
	}
    }

    printOut(\%prevGenes, $prevDNA);
    $prevDNA = $dnaSt;
    %prevGenes = %genes;
    $dnaSt = '';
    %genes = ();    
}

sub stripName
{
    my $in = shift;
    print "Error at $in !!! \n" if ($in eq '');
    # C095.3c.1
    my @name = split /\./,$in;
    $in = $name[0].$name[1] if (scalar @name == 3);
    $in =~ s/[a-z]$//;
    return $in;
}



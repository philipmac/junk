use strict;
use warnings;

## THIS OUTPUTS THE FOLLOWING ##
# Hits Reinke Expressed (>50) :647
# Hits Reinke Not Expressed (>50): 1730
# Intersection :0
# Whole PRG1 (hits plus no hits all expression levels):
 
# On I 2637 

# On II 3199 

# On III 2442 

# On IV 2914 

# On V 4554 

# On X 2544 
# Total Number of genes hit On IV  1715
# Total Number of genes hit NOT On IV  662
# Total Number of genes hit Not listed in DB at all  0
#  I : 156 
#  II : 133 
#  III : 121 
#  IV : 1715 
#   V : 173 
#  X : 79 
# Hits that were in the highly expressed genes Reinke Study:647

# Total Number of genes hit On IV + 169
# Total Number of genes hit NOT On IV + 478
#  I : 133 
#  II : 104 
#  III : 93 
#  IV : 169 
#   V : 104 
#  X : 44 
###

# # # # # # # # # # # # # # LOOK UP EVERYTHING IN PRG1 STUDY # 
my %wholePRG1Study;
open(IN, "/home/philip/other_data/prg1_alldata.clean.txt");
while (<IN>)
{
    # "T23D8.7" "T23D8.7"  "I" 1 9989791 9993503 1.1668253968254
    # C47B2.1" "C47B2.1" "fbxa-140" "I" 
    $_ =~ s/"//g;
    my (undef,$gene,undef,$csome, undef)= split /\s/,$_;
    $wholePRG1Study{$gene}=1;    
}
close IN;
# # # # # # # # # # # # # # 

my ($lowExpression_I,$lowExpression_II,$lowExpression_III,$lowExpression_IV,$lowExpression_V,$lowExpression_X);


my %genesNotOnIV;

open(IN, '/home/philip/other_data/I_Genes') or die $!;
my %genesOnI;
while (<IN>){
    chomp $_;
    next if $_ =~/^Seq/;
    $genesNotOnIV{$_}=1;
    $genesOnI{$_}=1;
    $lowExpression_I++ if $wholePRG1Study{$_};
}
close IN;

open(IN, '/home/philip/other_data/II_Genes') or die $!;
my %genesOnII;
while (<IN>){
    chomp $_;
    next if $_ =~/^Seq/;
    $genesOnII{$_}=1;
    $lowExpression_II++ if $wholePRG1Study{$_};
    $genesNotOnIV{$_}=1;
}
close IN;

open(IN, '/home/philip/other_data/III_Genes') or die $!;
my %genesOnIII;
while (<IN>){
    chomp $_;
    next if $_ =~/^Seq/;
    $genesOnIII{$_}=1;
    $lowExpression_III++ if $wholePRG1Study{$_};
    $genesNotOnIV{$_}=1;
}
close IN;

open(IN, '/home/philip/other_data/IV_Genes') or die $!;
my %genesOnIV;
while (<IN>){
    chomp $_;
    next if $_ =~/^Seq/;
    $lowExpression_IV++ if $wholePRG1Study{$_};
    $genesOnIV{$_}=1;
}
close IN;

open(IN, '/home/philip/other_data/V_Genes') or die $!;
my %genesOnV;
while (<IN>){
    chomp $_;
    next if $_ =~/^Seq/;
    $genesOnV{$_}=1;
    $lowExpression_V++ if $wholePRG1Study{$_};
    $genesNotOnIV{$_}=1;
}
close IN;

open(IN, '/home/philip/other_data/X_Genes') or die $!;
my %genesOnX;
while (<IN>){
    chomp $_;
    next if $_ =~/^Seq/;
    $genesOnX{$_}=1;
    $lowExpression_X++ if $wholePRG1Study{$_};
    $genesNotOnIV{$_}=1;
}
close IN;

open(IN, '/home/philip/other_data/bowtie/bin_offrate_1/cross_refd_outputs/all_hits_without_Reinke_expression') or die $!;
my %genesNotExpressed;


while (<IN>){
    chomp $_;
    my (undef, $gene, undef, undef) = split /\|/,$_;
    $genesNotExpressed{$gene}=1;
}
close IN;

open(IN, '/home/philip/other_data/bowtie/bin_offrate_1/cross_refd_outputs/all_hits_with_Reinke_expression') or die $!;
my %genesExpressed;
while (<IN>){
    chomp $_;
    my ($gene, undef, undef) = split /\t/,$_;
    $genesExpressed{$gene}=1;
}
close IN;


# ## ## ## ## ## ## ## ## ## ## ## ## ## #
# All hits #
my $totalOnIV;
my $totalNotOnIV;
my $genesNotListed =0;
my ($onI, $onII, $onIII, $onIV, $onV, $onX);

print 'Hits Reinke Expressed (>50) :';
print scalar keys %genesExpressed;

print "\nHits Reinke Not Expressed (>50): ";
print scalar keys %genesNotExpressed;
print "\n";

# print "Genes in PRG1 study, but <50: $lowExpression \n";

my $intersect = 0;
foreach (keys  %genesNotExpressed)
{
    $intersect++ if $genesExpressed{$_};
}
foreach (keys  %genesExpressed)
{
    $intersect++ if $genesNotExpressed{$_};
}

print "Intersection :$intersect\n";

print "Whole PRG1 (hits plus no hits all expression levels):\n 
On I $lowExpression_I \n
On II $lowExpression_II \n
On III $lowExpression_III \n
On IV $lowExpression_IV \n
On V $lowExpression_V \n
On X $lowExpression_X \n";

foreach (keys %genesNotExpressed, keys %genesExpressed){
    if($genesOnIV{$_}){ 
	$totalOnIV++ 
    }
    elsif ($genesNotOnIV{$_}) {
	$totalNotOnIV++
    }
    else {	
	print "cant get $_\n";
	$genesNotListed++;	# ie the guys that we just fail to find in the DB, should be zero idealy.
    }
    $onI++ if ($genesOnI{$_});
    $onII++ if ($genesOnII{$_});
    $onIII++ if ($genesOnIII{$_});
    $onIV++ if ($genesOnIV{$_});
    $onV++ if ($genesOnV{$_});
    $onX++ if ($genesOnX{$_});

}

print "Total Number of genes hit On IV  $totalOnIV\n";
print "Total Number of genes hit NOT On IV  $totalNotOnIV\n";
print "Total Number of genes hit Not listed in DB at all  $genesNotListed\n";

print " I : $onI \n II : $onII \n III : $onIII \n IV : $onIV \n  V : $onV \n X : $onX \n";


# ## ## ## ## ## ## ## ## ## ## ## ## ## #
# Reinke Expressed #
print "Hits that were in the highly expressed genes Reinke Study:";
print scalar keys %genesExpressed;
print "\n\n";

my ($onI_E, $onII_E, $onIII_E, $onIV_E, $onV_E, $onX_E);
my ($totalOnIV_E, $totalNotOnIV_E);
foreach (keys %genesExpressed){
    if($genesOnIV{$_}){ 
	$totalOnIV_E++ 
    }
    else {
	$totalNotOnIV_E++
    }

    $onI_E++ if ($genesOnI{$_});
    $onII_E++ if ($genesOnII{$_});
    $onIII_E++ if ($genesOnIII{$_});
    $onIV_E++ if ($genesOnIV{$_});
    $onV_E++ if ($genesOnV{$_});
    $onX_E++ if ($genesOnX{$_});
}

print "Total Number of genes hit On IV + $totalOnIV_E\n";
print "Total Number of genes hit NOT On IV + $totalNotOnIV_E\n";

print " I : $onI_E \n II : $onII_E \n III : $onIII_E \n IV : $onIV_E \n  V : $onV_E \n X : $onX_E \n";



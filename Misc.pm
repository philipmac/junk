package Misc;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(log_);

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
1;

#!/usr/bin/perl
use strict;
use warnings;

my @infiles = @ARGV;
print "file\tw non-venomous branches\tw venomous branches\n";
for my $file ( @infiles ) {
    open my $fh, '<', $file or die $!;
    while(<$fh>) {
        chomp;
        if ( /w \(dN\/dS\) for branches:  (\S+) (\S+)/ ) {
            my ( $nv, $v ) = ( $1, $2 );
            print $file, "\t", $nv, "\t", $v, "\n";
        }
    }
}
#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Bio::Phylo::Factory;
use Bio::Phylo::Cobra::TaxaMap;

my $TotalNchar = 0;

# get command line arguments
my ( $outgroup, $dir, $csv );
GetOptions(
    'dir=s' => \$dir,
    'csv=s' => \$csv,
    'outgroup=s' => \$outgroup,
);

# create taxa map
#my $map = Bio::Phylo::Cobra::TaxaMap->new($csv);

# concatenate tables
my %table = ( $outgroup => undef );
opendir my $dh, $dir or die $!;
while( my $entry = readdir $dh ) {
    if ( $entry =~ /\.dat/ ) {
        warn $entry;
        
        # going to read a simple key/value table
        open my $fh, '<', "$dir/$entry" or die $!;
		my $nchar;
		while(<$fh>) {
			chomp;
			my ( $key, $value ) = split /\t/, $_;
			$table{$key} .= $value;
			$nchar = length($value) if not defined $nchar;
			die "$entry - $nchar: $key $value" if $nchar != length($value);
		}
		$TotalNchar += $nchar;
    }
}

# create and populate matrix object
my $fac = Bio::Phylo::Factory->new;
my $matrix = $fac->create_matrix( '-type' => 'standard' );

# alphabetize, as taxa blocks also are
for my $code ( sort { $a cmp $b } keys %table ) {
    #my $binomial = $map->get_binomial_for_taxonID($row) || $outgroup;
    $matrix->insert(
    
    	# create and insert matrix row
        $fac->create_datum(
            '-type' => 'standard',
            '-name' => $code,
            '-char' => $code eq $outgroup ? '0' x $TotalNchar : $table{$code},
        )
    )
}

# sort rows alphabetically, as is done with taxa blocks
my @rows = sort { $a->get_name cmp $b->get_name } @{ $matrix->get_entities };
$matrix->clear;
$matrix->insert($_) for @rows; 

# print output
print "#NEXUS\n", $matrix->make_taxa->to_nexus, $matrix->to_nexus;

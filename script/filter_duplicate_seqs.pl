#!/usr/bin/perl
use strict;
use warnings;
use Bio::Phylo::IO 'parse_matrix';

my $matrix = parse_matrix(
	'-format' => 'fasta',
	'-type'   => 'dna',
	'-handle' => \*STDIN,
);

my %seen;
for my $row ( @{ $matrix->get_entities } ) {
	if ( my $name = $row->get_name ) {
		if ( not $seen{$name} ) {
			my $char = $row->get_char;
			print '>', $name, "\n", $char, "\n";
		}
		$seen{$name}++;
	}
}


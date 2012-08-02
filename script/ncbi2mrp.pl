#!/usr/bin/perl
use strict;
use Getopt::Long;
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse';
use Bio::Phylo::Cobra::TaxaMap;
use Bio::Phylo::Util::Logger;
use Bio::Phylo::Util::CONSTANT ':objecttypes';

# process command line arguments
my ( $infile, $format, $csv, $verbosity );
GetOptions(
    'infile=s' => \$infile,
    'format=s' => \$format,
    'csv=s'    => \$csv,
    'verbose+' => \$verbosity,
);

# instantiate logger
my $log = Bio::Phylo::Util::Logger->new( '-class' => 'main', '-level' => $verbosity );

# instantiate mapping
my $map = Bio::Phylo::Cobra::TaxaMap->new($csv);

# instantiate factory
my $fac = Bio::Phylo::Factory->new;

# here we parse as a forest, because forest objects
# can be turned in to MRP matrices
my ($forest) = @{
    parse(
        '-format' => $format,
        '-file'   => $infile,
        '-as_project' => 1,
    )->get_items(_FOREST_)
};
my $tree = $forest->first;

# here we map all binomials into species codes
my @binomials = $map->get_distinct_binomials;
for my $binomial ( @binomials ) {
	my $node;
	my $code = $map->get_code_for_binomial($binomial);
	
	# matches the entire name
	my $regex_full = qr/^\Q'$binomial'\E$/;
	
	# matches first parts (so could match subspecies)
	my $regex_partial = qr/^\Q'$binomial\E/;
	
	# these are just to shorten the subsequent invocation
	my %param  = ( '-value' => 'get_name' );
	my $method = 'get_by_regular_expression';
	
	# we have an exact match, simply rename to code
	if ( $node = $tree->$method( %param, '-match' => $regex_full )->[0] ) {
		$node->set_name( $code );
	}
	
	# we match a subspecies, which we need to give a sibling
	elsif ( $node = $tree->$method( %param, '-match' => $regex_partial )->[0] ) {
		my $sibling = $fac->create_node( '-name' => $code );
		$node->get_parent->set_child($sibling);
		$tree->insert($sibling);
	}
	else {
		$log->warn("NCBI taxonomy has no node for $binomial");
	}
}

# create MRP matrix
my $matrix = $forest->make_matrix;

# print MRP matrix as tab-delimited table
for my $row ( @{ $matrix->get_entities } ) {
    my $name = $row->get_name;
    my @char = $row->get_char;
    print $name, "\t", join('', @char), "\n";
}
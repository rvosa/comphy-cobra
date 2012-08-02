#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse';
use Bio::Phylo::Cobra::TaxaMap;
use Bio::Phylo::Util::Logger;
use Bio::Phylo::Util::CONSTANT ':objecttypes';

# process command line arguments
my ( $outgroup, $constr, $format, $ratchet, $csv, $verbosity );
GetOptions(
    'outgroup=s'   => \$outgroup,
    'constraint=s' => \$constr,
    'format=s'     => \$format,
    'ratchet=s'    => \$ratchet,
    'csv=s'        => \$csv,
    'verbose+'     => \$verbosity,
);

# instantiate helper objects
my $map = Bio::Phylo::Cobra::TaxaMap->new($csv);
my $log = Bio::Phylo::Util::Logger->new( '-class' => 'main', '-level' => $verbosity );
my $fac = Bio::Phylo::Factory->new;

# parse constraint tree (typically ncbi tree)
my ($tree) = @{ parse(
    '-format' => $format,
    '-file'   => $constr,
    '-as_project' => 1,
)->get_items(_TREE_) };

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

# remove those :4 branch lengths that NCBI adds for some reason
$tree->visit(sub{shift->set_branch_length(undef)});

# serialize back to newick
my $newick = $tree->to_newick;

# print command blocks
print <<"FOOTER";
BEGIN ASSUMPTIONS;
		OPTIONS DEFTYPE = ORD;
END;

BEGIN PAUP;
		OUTGROUP ${outgroup};
		CONSTRAINTS ncbi (BACKBONE) = ${newick}
		EXECUTE ${ratchet};
        QUIT;
END;
FOOTER

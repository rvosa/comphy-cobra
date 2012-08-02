#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Cobra::TaxaMap;

# process command line arguments
my ( $csv, $infile );
GetOptions(
	'csv=s'    => \$csv,
	'infile=s' => \$infile,
);

# instantiate map
my $map = Bio::Phylo::Cobra::TaxaMap->new($csv);

# parse tree
my $tree = parse_tree(
	'-format' => 'newick',
	'-file'   => $infile,
);

# traverse, label venomous branches
$tree->visit_depth_first(
	'-post' => sub {
		my $node = shift;
		if ( $node->is_terminal ) {
			my $name = $node->get_name;
			my $venom = $map->get_venom_for_phylip($name);
			$node->set_generic( 'venom' =>  $venom );
			$node->set_name( "${name}#${venom}" );
		}
		else {
			my @children = @{ $node->get_children };
			my @venomous = grep { $_->get_generic('venom') } @children;
			$node->set_generic( 'venom' => (@children == @venomous) );
			if ( $node->get_generic('venom') ) {
				$node->set_name('#1');
			}
			else {
				$node->set_name('#0');
			}
		}
	}
);

# print output
print 1, "\n", $tree->to_newick( '-nodelabels' => 1 );
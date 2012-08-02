#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::Phylo::Util::Logger;
use Bio::Phylo::Cobra::TaxaMap;
use Bio::Phylo::IO 'parse_tree';

# process command line arguments
my ( $infile, $csv, @skip, $format, $label );
my $verbose = 2;
GetOptions(
	'infile=s' => \$infile,
	'format=s' => \$format,
	'label=s'  => \$label,
	'csv=s'    => \$csv,
	'verbose+' => \$verbose,
	'skip=s'   => \@skip,
);
my %ignore = map { $_ => 1 } @skip;

my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbose,
	'-class' => 'main'
);

my $tree = parse_tree(
	'-format' => $format,
	'-file'   => $infile,
);

my $map = Bio::Phylo::Cobra::TaxaMap->new($csv);

my %tips_to_prune;
$tree->visit_depth_first(
	'-pre' => sub {
		my $node = shift;
		my $name = $node->get_name;
		
		# build up the path from root to tips
		if ( my $parent = $node->get_parent ) {
			my $path = ( $parent->get_generic('path') || 0 ) + ( $node->get_branch_length || 0 );
			$node->set_generic( 'path' => $path );
			$log->debug("path for $name: $path");
		}
		
		# start assembling distinct species sets for nodes
		if ( $node->is_terminal ) {
			my $method = 'get_code_for_' . $label;
			my $code = $map->$method( $node->get_name );
			$node->set_generic( 'species' => { $code => [ $node ] } );
			$log->debug("$code - $name");
		}
	},
	'-post' => sub {
		my $node = shift;		
		my @children = @{ $node->get_children };
		
		# node is internal
		if ( @children ) {
			$log->debug($node->get_name . ' is internal');
			
			my %species; # map of all subtended species codes and associated tips
			my @mononodes; # list of lists of monophyletic tips
			
			# iterate over children
			for my $child ( @children ) {
				my %child_species = %{ $child->get_generic('species') };
				my @codes = keys %child_species;
				
				# populate combined mapping
				for my $code ( @codes ) {
					$species{$code} = [] if not $species{$code};
					push @{ $species{$code} }, @{ $child_species{$code} };
				}
				
				# one key, multiple values
				if ( 1 == scalar @codes ) {
					my $code = shift @codes;
					if ( 1 < scalar @{ $child_species{$code} } && not $ignore{$code} ) {
						$log->info($child->get_name . ' is monophyletic: ' . Dumper(\%child_species));
						push @mononodes, [ @{ $child_species{$code} } ];
					}
				}
			}
			$node->set_generic( 'species' => \%species );
			
			# this node is not monophyletic, but some of its children are
			if ( scalar(keys %species) > 1 && scalar(@mononodes) > 0 ) {
				for my $set ( @mononodes ) {
					my @sorted = sort { $a->get_generic('path') <=> $b->get_generic('path') } @{ $set };
					my $nearest = shift @sorted;
					$tips_to_prune{$_->get_id} = $_ for @sorted;
				}
			}
		}
	}
);

$log->info("going to prune: ". $_->get_name) for values %tips_to_prune;
$tree->prune_tips([ values %tips_to_prune ]);
my %mapping;
$tree->visit(sub{
	my $node = shift;
	
	# start assembling distinct species sets for nodes
	if ( $node->is_terminal ) {
		my $name = $node->get_name;
		my $method = 'get_code_for_' . $label;
		my $code = $map->$method($name);
		$mapping{$code} = [] if not $mapping{$code};
		push @{ $mapping{$code} }, $name;
	}	
});
$log->info(Dumper(\%mapping));
print $tree->to_newick;
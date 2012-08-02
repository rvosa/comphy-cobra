use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Factory;
use Bio::Phylo::Cobra::TaxaMap;
use Bio::Phylo::IO qw'parse unparse';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Util::CONSTANT qw':objecttypes :namespaces';
use Data::Dumper;

# process command line arguments
my ( $infile, $format, $csv, $extension, $verbosity );
GetOptions(
    'stem=s'      => \$infile,
    'format=s'    => \$format,
    'csv=s'       => \$csv,
    'extension=s' => \$extension,
);

# instantiate helper objects
my $fac = Bio::Phylo::Factory->new;
my $log = Bio::Phylo::Util::Logger->new;
my $map = Bio::Phylo::Cobra::TaxaMap->new($csv);

# parse infile
my $project = parse(
    '-format'     => $format,
    '-file'       => "$infile.$extension",
    '-as_project' => 1
);

# get tree block from project
my ($forest) = @{ $project->get_items(_FOREST_) };

# resolve basal trichotomy, this is because PHYML writes NEWICK tree
# descriptions with three children at the root, which archeopteryx
# interprets as not fully resolved.
my ($tree) = @{ $forest->get_entities };
my $root = $tree->get_root;
my @children = @{ $root->get_children };
my $right1 = pop @children;
my $right2 = pop @children;
my $newroot = $fac->create_node;
$right1->set_parent($newroot);
$right2->set_parent($newroot);
$newroot->set_parent($root);
$tree->insert($newroot);

# make or fetch taxa block for trees block
my $taxa = $forest->make_taxa;
$project->insert($taxa);

# iterate over taxa
for my $taxon ( @{ $taxa->get_entities } ) {
    my $phylip = $taxon->get_name;
    
    my $code     = $map->get_code_for_phylip($phylip);
    my $binomial = $map->get_binomial_for_phylip($phylip);
    my $label    = $map->get_label_for_phylip($phylip);
            
	# attach scientific name and code as phyloxml annotations
	my %ns = ( 'pxml' => _NS_PHYLOXML_ );
	update_meta( $taxon, 'pxml:code' => $code, %ns );
	update_meta( $taxon, 'pxml:scientific_name' => $binomial, %ns );
	
	# use original sequence label as node name
	$taxon->get_nodes->[0]->set_name($label) if $label;
}

print unparse( '-format' => 'phyloxml', '-phylo' => $project );

# helper subroutine that attaches $predicate => $value to $object
# with namespace(s) %ns
sub update_meta {
    my ( $object, $predicate, $value, %ns ) = @_;
    my ( $meta ) = @{ $object->get_meta($predicate) };
    if ( $meta ) {
        $meta->set_triple( $predicate => $value );
    }
    else {
        $object->add_meta(
            $fac->create_meta(
                '-namespaces' => \%ns,
                '-triple' => { $predicate => $value },
            )
        );
    }
}
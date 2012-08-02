#!/usr/bin/perl
use strict;
use Getopt::Long;
use Bio::Phylo::Factory;
use Bio::Phylo::IO qw'parse parse_matrix';
use Bio::Phylo::Util::CONSTANT ':objecttypes';

# instantiate factory
my $fac = Bio::Phylo::Factory->new;

# process command line arguments
my ( $treefile, $treeformat, $datafile, $dataformat, $datatype, $labels );
GetOptions(
    'treefile=s'   => \$treefile,
    'treeformat=s' => \$treeformat,
    'datafile=s'   => \$datafile,
    'dataformat=s' => \$dataformat,
    'datatype=s'   => \$datatype,
    'labels+'      => \$labels,
);

# read trees from tree file
my ($forest) = @{
    parse(
        '-format' => $treeformat,
        '-file'   => $treefile,
        '-as_project' => 1,
    )->get_items(_FOREST_)
};

# optionally apply node labels for hyphy branchsiterel
my %tips;
if ( $labels ) {
    $forest->visit(sub{
        my $tree = shift;
        my $i = 0;
        $tree->visit(sub{
            my $node = shift;
            if ( $node->is_internal ) {
                $node->set_name( 'node' . ++$i );
            }
            else {
                $tips{$node->get_name} = 1;
            }
        });
    });
}

# read data from data file
my $matrix = parse_matrix(
    '-format' => $dataformat,
    '-file'   => $datafile,
    '-type'   => $datatype,
);

# prune sequences not in tree
my @delete;
$matrix->visit(sub{
    my $row = shift;
    push @delete, $row if not $tips{$row->get_name};
});
$matrix->delete($_) for @delete;

# reconcile taxa
my $proj = $fac->create_project;
my $treetaxa = $forest->make_taxa;
my $datataxa = $matrix->make_taxa;
my $merged = $treetaxa->merge_by_name($datataxa);
$forest->set_taxa($merged);
$matrix->set_taxa($merged);
my @sorted = sort { $a->get_name cmp $b->get_name } @{ $matrix->get_entities };
$matrix->clear;
$matrix->insert($_) for @sorted;

# print result
print $proj->insert($merged,$forest,$matrix)->to_nexus( -nodelabels => 1 );

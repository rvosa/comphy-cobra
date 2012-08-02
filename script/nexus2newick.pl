#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO qw'parse';
use Bio::Phylo::Util::CONSTANT ':objecttypes';

# process command line arguments
my ( $infile );
GetOptions(
    'infile=s'  => \$infile,
);

# fetch first tree
my ($tree) = @{
    parse(
        '-format' => 'nexus',
        '-file'   => $infile,
        '-as_project' => 1,
    )->get_items(_TREE_)
};

# print output
print $tree->to_newick;
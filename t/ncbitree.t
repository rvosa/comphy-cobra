#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';
use Bio::Phylo::IO 'parse';
use Bio::Phylo::Cobra::TaxaMap;

my $script = $0;
my $csv = $script;
my $phyliptree = $script;
$csv =~ s|ncbitree\.t$|../data/excel/taxa.csv|;
$phyliptree =~ s|ncbitree\.t$|../data/phyliptree.phy|;

my $map = Bio::Phylo::Cobra::TaxaMap->new($csv);
my $tree = parse( '-format' => 'newick', '-file' => $phyliptree )->first;
my %seen_in_csv = map { $_ => 1 } $map->get_distinct_codes;
my %seen_in_tree;
for my $tip ( @{ $tree->get_terminals } ) {
    my $name = $tip->get_name;
    $name =~ s/'//g;
    $seen_in_tree{$map->get_code_for_binomial($name)} = 1;
}

# test to see that all csv taxa are in the ncbi common tree
for my $name ( keys %seen_in_csv ) {
    ok( $seen_in_tree{$name}, "$name from csv is in ncbi common tree" );
}

# test to see that all tips in the common tree are in the csv
for my $name ( keys %seen_in_tree ) {
    ok( $seen_in_csv{$name}, "$name from common tree is in csv" );
}
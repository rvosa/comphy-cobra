#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';
use Bio::Phylo::Cobra::TaxaMap;

my $csv = $0;
$csv =~ s|speciescodes\.t$|../data/excel/taxa.csv|;
my $map = Bio::Phylo::Cobra::TaxaMap->new($csv);

my @names = $map->get_distinct_binomials;
my %variants;
for my $name ( @names ) {
    my @parts = split /\s/, $name;
    my $binomial = $parts[0] . ' ' . $parts[1];
    $variants{$binomial} = [] unless $variants{$binomial};
    push @{ $variants{$binomial} }, $name;
}

for my $taxon ( keys %variants ) {
    my @across_variants;
    for my $variant ( @{ $variants{$taxon} } ) {
        my @codes = $map->get_code_for_binomial($variant);
        my @uniq = keys %{ { map { $_ => 1 } @codes } };
        
        # test that there is only one code per scientific name
        ok( scalar(@uniq) == 1, "@uniq" );
        push @across_variants, @uniq;
    }
    my @uniq = keys %{ { map { $_ => 1 } @across_variants } };
    
    # test that there is only one code for all subspecies
    ok( scalar(@uniq) == 1, "@across_variants" );
    
    my %seen = map { $_ => 1 } @{ $variants{$taxon} };
    for my $code ( @across_variants ) {
        my @names = $map->get_binomial_for_code($code);
        for my $name ( @names ) {
            
            # test that the code is not used for multiple binomials
            ok( $seen{$name}, "$name" );
        }
    }
}


my @phylips = $map->get_all_phylips;
my @uniq = keys %{ { map { $_ => 1 } @phylips } };
    
# test that phylip names are unique per gene
ok(scalar(@uniq) == scalar(@phylips), "check phylip names are GUIDs");
    
for my $phylip ( @phylips ) {
    my $code = $map->get_code_for_phylip($phylip);
    
    # test that phylip names are codes with optional numerical suffix
    ok( $phylip =~ /^$code\d*$/, "$phylip => $code" );
}
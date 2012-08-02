#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Cobra::TaxaMap;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::IO qw'parse_matrix unparse';

# process command line arguments
my ( $csv, @ignore, $verbose );
my $threshold = 0.9;
GetOptions(
    'csv=s'       => \$csv,
    'ignore=s'    => \@ignore,
    'threshold=s' => \$threshold,
    'verbose+'    => \$verbose,
);
my %ignore = map { $_ => 1 } @ignore;

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new( '-level' => $verbose, '-class' => 'main' );
my $map = Bio::Phylo::Cobra::TaxaMap->new( $csv );

# parse matrix from stdin
my $matrix = parse_matrix(
    '-format' => 'fasta',
    '-handle' => \*STDIN,
    '-type'   => 'dna'
);

# cluster rows by species code
my %rows_by_code;
$matrix->visit(
    sub {
        my $row   = shift;
        my $name  = $row->get_name;
        my $label = $map->parse_label($name);
        my $code  = $map->code($label);
        $rows_by_code{$code} = [] unless $rows_by_code{$code};
        push @{ $rows_by_code{$code} }, $row;
        $row->set_generic( 'fasta_def_line' => $row->get_name );
    }
);

# filter short rows
CLUSTER: for my $code ( keys %rows_by_code ) {
    
    # this is so we don't filter cobra or anolis, which we think are
    # true paralogs
    if ( $ignore{$code} ) {
        $log->debug("ignoring cluster $code");
        next CLUSTER;
    }
    $log->debug("filtering cluster $code");
    
    # compute longest seq length
    my $longest = 0;
    for my $row ( @{ $rows_by_code{$code} } ) {
        my $seq = $row->get_unaligned_char;
        $longest = length($seq) if length($seq) > $longest;
    }
    $log->debug("longest seq is $longest bp");
    
    # now filter out those in the cluster whose length is less than
    # the threshold relative to the longest
    for my $row ( @{ $rows_by_code{$code} } ) {
        my $seq = $row->get_unaligned_char;
        my $length = length($seq);
        if ( $length / $longest < $threshold ) {
            $log->debug("deleting row " . $row->get_name);
            $matrix->delete($row);
        }
    }
}

# done
print unparse( '-format' => 'fasta', '-phylo' => $matrix );

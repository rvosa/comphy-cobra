#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::IO qw'parse_matrix unparse';

# process command line arguments
my $verbosity;
my $cutoff = 0.5;
GetOptions(
    'cutoff=s' => \$cutoff,
    'verbose+' => \$verbosity,
);

# instantiate logger
my $log = Bio::Phylo::Util::Logger->new( '-class' => 'main', '-level' => $verbosity );

# parse input matrix
my $matrix = parse_matrix(
    '-format' => 'fasta',
    '-handle' => \*STDIN,
    '-type'   => 'dna',
);

# start iterating over the columns
my $nchar = $matrix->get_nchar;
my $ntax  = $matrix->get_ntax;
my $raw   = $matrix->get_raw;

# this will hold zero-based indices of columns to remove
my @indices;
for ( my $i = 1; $i <= $nchar; $i += 3 ) {
    $log->debug("checking site $i");
    
    # counts rows with missing data for first codon position
    my $count = 0;
    for my $j ( 0 .. ( $ntax - 1 ) ) {
        
        # only checks first codon position
        my $state = $raw->[$j]->[$i];
        if ( $state eq '?' or $state eq '-' ) {
            my $name = $raw->[$j]->[0];
            $log->debug("$name has missing/gap state $state at site $i");
            $count++;
        }
    }
    
    # compute fraction missing/gap
    my $fraction = $count / $ntax;
    $log->debug("fraction of missing at site $i is $fraction");    
    if ( $fraction > $cutoff ) {
        
        # watch out: raw matrix's first field is row name
        push @indices, ( $i - 1 );
        push @indices, $i;
        push @indices, ( $i + 1 );
    }
}

# make taxa so that the cloning works
my $taxa = $matrix->make_taxa;
$matrix->set_taxa($taxa);

# do the deletion
my $pruned = $matrix->prune_chars(\@indices);

# rename fasta_def_line or we get two >> at beginning
$pruned->visit(
    sub {
        my $row = shift;
        $row->set_generic( 'fasta_def_line' => $row->get_name );
    }
);

# return output
print unparse( '-format' => 'fasta', '-phylo' => $pruned );


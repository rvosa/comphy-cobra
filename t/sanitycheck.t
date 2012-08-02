#!/usr/bin/perl
use strict;
use Test::More 'no_plan';
use Bio::Phylo::IO 'parse';
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Data::Dumper;

my $dir = $0;
$dir =~ s|sanitycheck\.t$|../data/raw/|;

warn "Going to read data from $dir\n";

opendir my $dirhandle, $dir or die $!;
while( my $entry = readdir $dirhandle ) {
    
    # for every *.fas file there must exist a *.nex file
    if ( $entry =~ m/(.+)\.fas$/ ) {
        my $stem = $1;
        my $treefile = $dir . $stem . '.tre';
        my $seqfile  = $dir . $stem . '.fas';
        ok( -e $treefile && -e $seqfile, "$treefile and $seqfile exist");
        
        # now parse the sequence file
        my ($seqs) = @{ parse(
            '-format' => 'fasta',
            '-type'   => 'dna',
            '-file'   => $seqfile,
            '-as_project' => 1,
        )->get_items(_MATRIX_) };
        
        # now parse the tree file
        my @trees = @{ parse(
            '-format' => 'nexus',
            '-file'   => $treefile,
            '-as_project' => 1,
        )->get_items(_TREE_) };
        
        # test if numbers match up
        my $seqntax = $seqs->get_ntax;
        for my $tree ( @trees ) {
            ok( $tree->get_ntax == $seqntax, "numbers match in $stem" );
        }
        
        # populate %seqnames and %tipnames hashes with first 3 parts of names
        my $pattern = qr/^([^_]+_[^_]+_[^_]+)/;
        my ( %seqnames, %tipnames );
        for my $seq ( @{ $seqs->get_entities } ) {
            if ( $seq->get_name =~ $pattern ) {
                my $name = $1;
                $seqnames{$name}++;
                ok( $seqnames{$name} == 1, "$name not seen yet in $seqfile" );
            }
            else {
                ok( undef, $seq->get_name );
            }
        }

        for my $tree ( @trees ) {
            for my $tip ( @{ $tree->get_terminals } ) {
                if ( $tip->get_name =~ $pattern ) {
                    my $name = $1;
                    $tipnames{$name}++;
                    ok( $tipnames{$name} <= 2, "$name seen <= 2 times in $treefile" );
                }
                else {
                    ok( undef, $tip->get_name );
                }            
            }
        }
        
        # test if the sequence names occur in the tree file
        for my $seqname ( keys %seqnames ) {
            ok( $tipnames{$seqname}, "seq $seqname occurs in $treefile" );
        }
        for my $tipname ( keys %tipnames ) {
            ok( $seqnames{$tipname}, "tip $tipname occurs in $seqfile" );
        }        
        
    }
}



#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Test::More 'no_plan';
use Bio::DB::GenBank;
use Bio::Phylo::IO 'parse';
use Bio::Tools::CodonTable;
use Bio::Phylo::Cobra::TaxaMap;
use Bio::Phylo::Util::CONSTANT ':objecttypes';

# compute paths
my $dir = $0;
$dir =~ s|proteintest\.t|../data/raw|;
my $csv = $dir;
$csv =~ s|raw$|excel/taxa.csv|;

# instantiate objects
my $map   = Bio::Phylo::Cobra::TaxaMap->new($csv);
my $table = Bio::Tools::CodonTable->new;
my $dbgb  = Bio::DB::GenBank->new;

# open dir handle to fasta file dir
opendir my $dh, $dir or die $!;
while( my $entry = readdir $dh ) {
    if ( $entry =~ /\.fas/ ) {
        
        # parse fasta file
        my ($matrix) = @{
            parse(
                '-format'     => 'fasta',
                '-file'       => "$dir/$entry",
                '-as_project' => 1,
                '-type'       => 'dna',
            )->get_items(_MATRIX_)
        };
        
        # create amino acid translation
        for my $row ( @{ $matrix->get_entities } ) {
            my @char = $row->get_char;
            my $translation = '';
            for ( my $i = 0; $i <= $#char; $i += 3 ) {
                my $codon = $char[$i] . $char[$i+1] . $char[$i+2];
                my $aa = $table->translate($codon);
                $translation .= $aa if $aa && $aa ne 'X' && $aa ne '-';
            }
            
            # get protein from genbank
            my $label = $map->parse_label($row->get_name);
            if ( my $gi = $map->gi($label) ) {
                eval {
                    my $seq = $dbgb->get_Seq_by_gi($gi);
                    
                    # get the coding sequence	
                    for my $feat ( $seq->get_SeqFeatures ) {
                        if ( $feat->primary_tag eq 'CDS' ) {
                            my ($feat_seq) = @{ $feat->{'_gsf_tag_hash'}->{'translation'} };
                            ok( $feat_seq =~ /\Q$translation\E/ || $translation =~ /\Q$feat_seq\E/, "$entry => $label" );
                            if ( $feat_seq !~ /\Q$translation\E/ && $translation !~ /\Q$feat_seq\E/ ) {
                                warn "Genbank: $feat_seq\nFasta:   $translation\n";    
                            }
                        }
                    }
                };
                warn $@ if $@;
            }
        }
    }
}
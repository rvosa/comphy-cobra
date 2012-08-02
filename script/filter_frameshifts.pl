#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::DB::GenBank;
use Bio::Phylo::Factory;
use Bio::Tools::CodonTable;
use Bio::Phylo::Cobra::TaxaMap;
use Bio::Phylo::IO qw'parse unparse';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Util::CONSTANT qw':objecttypes :namespaces';

# process command line arguments
my ( $csv );
GetOptions(
    'csv=s' => \$csv,
);

# instantiate helper objects
my $dbgb  = Bio::DB::GenBank->new;
my $table = Bio::Tools::CodonTable->new;
my $map   = Bio::Phylo::Cobra::TaxaMap->new($csv);
my $log   = Bio::Phylo::Util::Logger->new(
    '-level' => INFO,
    '-class' => 'main'
);

# parse matrix from fasta alignment
my ($matrix) = @{
    parse(
        '-format' => 'fasta',
        '-handle' => \*STDIN,
        '-type'   => 'dna',
        '-as_project' => 1,
    )->get_items(_MATRIX_)
};

# filter incorrect back translations
my @incorrect;
ROW: for my $row ( @{ $matrix->get_entities } ) {
    my @char = $row->get_char;
    
    # compute protein translations
    my $translation = '';
    for ( my $i = 0; $i <= $#char; $i += 3 ) {
        my $codon = $char[$i] . $char[$i+1] . $char[$i+2];
        my $aa = $table->translate($codon);
        $translation .= $aa if $aa && $aa ne 'X' && $aa ne '-';
    }
    $log->debug("computed protein translation");
    
    # get protein from genbank
    my $name = $row->get_name;
    $name =~ s/^>//;
    $row->set_generic( 'fasta_def_line' => $name );
    my $label = $map->parse_label($name);
    if ( my $gi = $map->gi($label) ) {
        $log->debug("going to verify $gi");
        eval {
            my $seq = $dbgb->get_Seq_by_gi($gi);
            
            # get the coding sequence	
            for my $feat ( $seq->get_SeqFeatures ) {
                if ( $feat->primary_tag eq 'CDS' ) {                    
                    my ($feat_seq) = @{ $feat->{'_gsf_tag_hash'}->{'translation'} };
                    $feat_seq =~ s/X//g;
                    if ( $feat_seq !~ /\Q$translation\E/ && $translation !~ /\Q$feat_seq\E/ ) {
                        $log->warn( "***Possible frame shift:\nGenbank: $feat_seq\nFasta:   $translation\n" );
                        push @incorrect, $row;
                    }
                    else {
                        $log->info("no frame shifts in alignment of $gi");
                        next ROW;
                    }
                }
            }
        };
    }    
}
for my $row ( @incorrect ) {
    $matrix->delete($row);
}

print unparse( '-format' => 'fasta', '-phylo' => $matrix );
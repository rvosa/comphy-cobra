#!/usr/bin/perl
use strict;
use Getopt::Long;
use Bio::Phylo::IO 'parse';
use Bio::Phylo::Cobra::TaxaMap;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Data::Dumper;

# process command line arguments
my ( $infile, $format, $csv, $verbose );
GetOptions(
    'infile=s' => \$infile,
    'format=s' => \$format,
    'csv=s'    => \$csv,
	'verbose+' => \$verbose,
);

my $log = Bio::Phylo::Util::Logger->new( '-level' => $verbose, '-class' => 'main' );

# create seen hash for phyloxml codes
my $map = Bio::Phylo::Cobra::TaxaMap->new($csv);
my %seen = map { $_ => 1 } $map->get_all_codes;
$log->debug(Dumper(\%seen));

# parse tree block from input file
my ($forest) = @{
    parse(
        '-format' => $format,
        '-file'   => $infile,
        '-as_project' => 1,
    )->get_items(_FOREST_)
};

# for each tip in each tree, fetch its taxon object and from that get the
# skos:*match annotations, which may have an NCBI taxon id. if it does,
# use it to find the phyloxml code and use that for tip and taxon name
for my $tree ( @{ $forest->get_entities } ) {
	$log->debug("tree: " . $tree->get_name() . " in file $infile");
	
    for my $tip ( @{ $tree->get_terminals } ) {
        my $taxon = $tip->get_taxon;
        my @predicates = qw(skos:closeMatch skos:exactMatch);
        if ( $taxon) {
			$log->debug("taxon: " . $taxon->get_name);
			
			# loop over metadata annotations that match @predicates
			META: for my $meta ( @{ $taxon->get_meta(@predicates) } ) {
				my $obj = $meta->get_object;
				$log->debug($obj);
				
				# this matches the ncbi taxonomy				
				if ( $obj =~ m|http://purl.uniprot.org/taxonomy/(\d+)| ) {
					my $id = $1;
					$log->debug($id);
					my $code = $map->get_code_for_taxonID($id);
					$tip->set_name($code);
					$taxon->set_name($code);
					last META;
				}
			}
        }
        else {
        	$log->warn("no taxon for tip " . $tip->get_xml_id . " in $infile");
        }
    }
}

# create mrp matrix
my $matrix = $forest->make_matrix;

# create a simple hash keyed on ncbi taxon ids, with values the character
# state sequences. only keep those key value pairs that are seen in taxa.csv
my %simple;
my $nchar;
for my $row ( @{ $matrix->get_entities } ) {
    my $name = $row->get_name;
	$log->debug('row: '.$name);
	
    if ( $seen{$name} ) {		
        my @char = $row->get_char;
		$log->debug($name."\t".join('',@char));
        $simple{$name} = \@char;
        $nchar = scalar @char;
    }
}

# only keep phylogenetically informative columns. these mrp matrices *can*
# have uninformative columns because we've pruned rows. also, they now may
# have duplicate 'site patterns', which we also prune
my %informative = map { $_ => [] } keys %simple;
my @names = keys %simple;
my %pattern;
for my $i ( 0 .. ( $nchar - 1 ) ) {
    my ( %char, @char, $pattern );    
    for my $name ( @names ) {
        $char{$simple{$name}->[$i]}++;
        push @char, $simple{$name}->[$i];
        $pattern .= $simple{$name}->[$i];
    }
    if ( scalar(keys(%char)) > 1 && ! $pattern{$pattern} ) {
        for my $j ( 0 .. $#names ) {
            push @{ $informative{$names[$j]} }, $char[$j];
        }
    }
    $pattern{$pattern}++;
}

# calculate pruned, informative, distinct nchar
my $informnchar;
ROW: for my $row ( keys %informative ) {
	if ( not defined $informnchar ) {
		$informnchar = scalar @{ $informative{$row} };
		next ROW;
	}
	die $informnchar if $informnchar != scalar @{ $informative{$row} };
}

# print as simple key/value table
for my $row ( keys %seen ) {
	if ( $informative{$row} ) {
		my $seq = join '', @{ $informative{$row} };
		$seq .=  '?' x ( $informnchar - length($seq) );
    	print $row, "\t", $seq, "\n";
    }
    else {
    	print $row, "\t", ( '?' x $informnchar ), "\n";
    }
}

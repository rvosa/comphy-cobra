#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger;
use Bio::Phylo::Cobra::TaxaMap;
use Bio::Phylo::IO qw'parse unparse';
use Bio::Phylo::Util::CONSTANT ':objecttypes';

# process command line arguments
my ( $infile, $serializer, $mrbayesblock, $csv, $datatype, $verbosity );
GetOptions(
	'infile=s'       => \$infile,
	'serializer=s'   => \$serializer,
	'mrbayesblock=s' => \$mrbayesblock,
	'datatype=s'     => \$datatype,
	'csv=s'          => \$csv,
	'verbose+'       => \$verbosity,
);

# instantiate logger
my $log = Bio::Phylo::Util::Logger->new( '-class' => 'main', '-level' => $verbosity );

# instantiate taxa map
my $map = Bio::Phylo::Cobra::TaxaMap->new($csv);

# parse data
my ($matrix) = @{parse(
	'-file'       => $infile,
	'-format'     => $serializer,
	'-type'       => $datatype,
	'-as_project' => 1,
)->get_items(_MATRIX_)};

# emit logging msg about alignment parsing
$log->info("parsed $serializer data of type $datatype from file $infile");

# rename tips
$matrix->visit(sub{
	my $row = shift;
	my $name = $row->get_name;
	
	# do the lookup and rename
	if ( my $phylip = $map->phylip($row->get_name) ) {
		$row->set_name($phylip);
		$log->info("renamed row $name to short name $phylip");
	}
	else {
		$log->error("couldn't find short phylip name for row $name");
	}
});

# produce output
print "#NEXUS\n";
print $matrix->to_nexus( '-data_block' => 1 );

# print mrbayes block, if any
if ( $mrbayesblock ) {
	if ( -e $mrbayesblock ) {
		$log->info("going to copy mrbayes block from file $mrbayesblock");
		open my $fh, '<', $mrbayesblock or die $!;
		while(<$fh>) {
			print $_;
		}
	}
	else {
		$log->info("mrbayes block provided as command line argument");
		print $mrbayesblock;
	}
}
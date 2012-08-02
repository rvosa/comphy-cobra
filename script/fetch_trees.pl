#!/usr/bin/perl
use strict;
use warnings;
use Bio::Phylo::Cobra::TaxaMap;
use Bio::Phylo::PhyloWS::Client;
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse';
use Bio::Phylo::Util::Logger;
use Getopt::Long;

# process command line arguments
my ( $dir, $file, $id, $verbosity );
GetOptions(
	'dir=s' => \$dir, # sourcetrees
	'csv=s' => \$file, # taxa.csv
	'id=i'  => \$id, # optional: use single NCBI identifier
	'verbose+' => \$verbosity,
);

# instantiate logger
my $log = Bio::Phylo::Util::Logger->new( '-class' => 'main', '-level' => $verbosity );

# read csv file, get distinct NCBI taxon IDs
my $map = Bio::Phylo::Cobra::TaxaMap->new($file);
my @ids = $map->get_distinct_taxonIDs;
@ids = ( $id ) if $id;
$log->debug("going to fetch trees for these NCBI taxon IDs: @ids");

# instantiate client
my $fac = Bio::Phylo::Factory->new;
my $client = $fac->create_client(
    '-base_uri'  => 'http://treebase.org/treebase-web/phylows/',
    '-authority' => 'TB2',
);

# fetch all trees for all identifiers
my %seen;
my @notseen;
for my $id ( @ids ) {
	eval {
		my $desc = $client->get_query_result(
			'-query'        => "tb.identifier.ncbi=$id",
			'-section'      => 'taxon',
			'-recordSchema' => 'tree',
		);
		
		# XXX the treebase guid is the purl
		for my $res ( @{ $desc->get_entities } ) {
			my $guid = $res->get_guid;
			
			# write the result to a file based on the GUID.
			my $outfile = $guid;
			$outfile =~ s|^.+:(.+)$|$dir/$1.xml|;				
			
			# don't download twice
			if ( not -e $outfile ) {
				
				# fetch data
				my $url  = $guid . '?format=nexml';
				my $proj = parse( 
					'-format'     => 'nexml', 
					'-url'        => $url,
					'-as_project' => 1,
				);						
				open my $fh, '>', $outfile or die "Can't open $outfile: $!";
				print $fh $proj->to_xml, "\n";
				$log->info("taxon $id written to $outfile");
			}
			else {
				$log->info("already downloaded $outfile");
			}
			
			# counting seen files so we can count overlap
			$seen{$guid}++;
		}
		
		# zero hits
		push @notseen, $id unless scalar( @{ $desc->get_entities } );
	};
	if ( $@ ) {
		$log->warn( "problem with taxon $id: $@" );
	}
}

# reporting back
my $ntax = scalar(@ids);
my $filecount = scalar(keys(%seen));
my $average = $ntax / $filecount;
print "$average overlapping taxa in $filecount files\n";
print "taxa with 0 hits:\n";
print join "\n", @notseen;

__DATA__
0.414860681114551 overlapping taxa in 323 files
taxa with 0 hits:
111304
186611
196418
310520
33626
338837
338838
51750
537493
61970
672774
865857
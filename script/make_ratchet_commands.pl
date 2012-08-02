#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_matrix';
use Bio::Phylo::Util::Logger;

# process command line arguments
my ( $setup, $matrix, $reps, $pct, $verbosity );
GetOptions(
    'setup=s'  => \$setup,
    'matrix=s' => \$matrix,
    'reps=i'   => \$reps,
    'pct=i'    => \$pct,
    'verbose+' => \$verbosity,
);

# instantiate logger
my $log = Bio::Phylo::Util::Logger->new( '-class' => 'main', '-level' => $verbosity );
$log->info("setup ... $setup");
$log->info("matrix .. $matrix");
$log->info("reps .... $reps");
$log->info("pct ..... $pct");

# parse matrix for nchar
my $nchar = parse_matrix('-format' => 'nexus', '-file' => $matrix)->get_nchar;
$log->info("parsed matrix $matrix, nchar is $nchar");

# will only read these lines from setup.nex
my %commands = (
    'startcmd' => [],
    'paupcmd'  => [],
    'normcmd'  => [],
    'stopcmd'  => [],
);

# read setup.nex
{
    $log->info("going to read ratchet template commands from $setup");
    open my $fh, '<', $setup or die $!;
    while(<$fh>) {
        chomp;
        if ( /(startcmd|paupcmd|normcmd|stopcmd)\s*"(.+)"/ ) {
            my ( $prefix, $line ) = ( $1, $2 );
            if ( exists $commands{$prefix} ) {
                push @{ $commands{$prefix} }, $line;
                $log->debug("adding $prefix - $line");
            }
        }
        else {
            $log->debug("skipping line $_");
        }
    }
}

# print header
print << 'HEADER';
#nexus
begin paup;
[**** starting commands ****]
HEADER

# print starting commands
print_commands('startcmd');
$log->info("printed file header and start commands");

# print iterations
for my $i ( 1 .. $reps ) {
    $log->debug("iteration: $i");
    print << "ITERATION";
[!
****************************
*** starting iteration $i ***
****************************
]
[*** reweighting step for iteration $i ***]
ITERATION

    # assign weights by sampling with replacement
    my %weight;
    for my $j ( 1 .. int( $nchar * ($pct/100) ) ) {
        my $char = int(rand($nchar+1));
        $weight{$char}++;
    }
    
    # bin by weight
    my %binned;
    for my $char ( 1 .. $nchar ) {
        my $weight = exists $weight{$char} ? $weight{$char} + 1 : 1; 
        $binned{$weight} = [] if not exists $binned{$weight};
        push @{ $binned{$weight} }, $char;
    }
    
    # get all weights
    my @weights = sort { $a <=> $b } keys %binned;
    
    # print weighting statements
    for my $weight ( @weights ) {
        print "weights ${weight}: ", join(' ',@{ $binned{$weight} }), ";\n";
    }
    
    # print search commands while re-weighted
    print_commands('paupcmd');
    
    # print unweighted commands
    print << "UNWEIGHTED";
[*** returning to original weights for iteration $i ***]
weights 1: 1-${nchar};
UNWEIGHTED

print_commands('paupcmd');
print_commands('normcmd');
}
$log->info("printed iterations");

# print stopping commands
print "[**** stopping commands ****]\n";
print_commands('stopcmd');
print "end;\n";
$log->info("printed stopping commands");

sub print_commands {
    my $key = shift;
    print "$_;\n" for @{ $commands{$key} };
}
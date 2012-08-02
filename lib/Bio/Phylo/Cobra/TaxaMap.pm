package Bio::Phylo::Cobra::TaxaMap;
use strict;
use warnings;
use Bio::Phylo::Util::Exceptions 'throw';
our $AUTOLOAD;

my %keys = (
	'label'    => 0,
	'binomial' => 1,
	'taxonID'  => 2,
	'code'     => 3,
	'gi'       => 4,
	'phylip'   => 5,
	'gene'     => 6,
	'venom'    => 7,
);

=pod

=head1 GETTER/SETTERS

=over

=item $map->binomial($label, [$binomial])

=item $map->taxonID($label, [$taxonID])

=item $map->gi($label, [$gi])

=item $map->code($label, [$code])

=item $map->phylip($label, [$phylip])

=item $map->gene($label, [$gene])

=item $map->venom($label, [$venom])

=back

=head1 QUERIES

In the following, x and y are both one of label, 
binomial, taxonID, gi, code, phylip, gene, venom.

=over

=item $map->get_x_for_y()

In scalar context, this returns the first result. In list context, all matching
x will be returned in case of a one-to-many cardinality between y and x.

=item $map->get_all_xs()

=item $map->get_distinct_xs()

=back

=head1 METHODS

=over

=item to_csv()

=item as_2d()

=back

=cut


sub new {
	my $package = shift;
	my $file = shift;
	my %self;
	open my $fh, '<', $file or die $!;
	while(<$fh>) {
		chomp;
		my @fields = split /,/, $_;
		my $key = shift @fields;
		$self{$key} = \@fields;
	}
	return bless \%self, $package;
}

sub parse_label {
	my ( $self, $label ) = @_;
	if ( $label =~ /^([^_]+_[^_]+_[^_]+)/ ) {
		my $parsed = $1;
		return $parsed;
	}
	throw 'BadArgs' => "Bad label '$label', expecting >= 3 underscore-separated words"
}

sub as_2d {
	my $self = shift;
	my @result;
	for my $key ( sort { $a cmp $b } keys %{ $self } ) {
		push @result, [ $key, @{ $self->{$key} } ];
	}
	return @result;
}

sub to_csv {
	my $self = shift;
	my $string = '';
	for my $row ( $self->as_2d ) {
		no warnings 'uninitialized'; # there *can* be empty fields, it's ok
		$string .= join ',', @{ $row };
		$string .= "\n";	
	}
	return $string;
}

sub AUTOLOAD {
	my $method = $AUTOLOAD;
	$method =~ s/.*://;
	if ( exists $keys{$method} ) {
		my $self = shift;
		my $label = shift;
		my $value = shift;
		if ( defined $value ) {
			# when doing a hash-based lookup, the indices
			# are off by one because $label isn't part of
			# the value array
			$self->{$label}->[$keys{$method} - 1] = $value;
		}
		return $self->{$label}->[$keys{$method} - 1];
	}
	elsif ( $method =~ m/get_([^_]+)_for_([^_]+)/ ) {
		my ( $wanted, $key ) = ( $1, $2 );
		my $self = shift;
		my $value = shift;
		my @result;
		for my $row ( $self->as_2d ) {
			if ( $row->[$keys{$key}] && $row->[$keys{$key}] eq $value ) {
				push @result, $row->[$keys{$wanted}];
			}
		}
		return wantarray ? @result : $result[0];
	}
	elsif ( $method =~ m/get_all_([^_]+)s/ ) {
		my $field = $1;
		my $self = shift;
		my @result;
		for my $row ( $self->as_2d ) {
			push @result, $row->[$keys{$field}];
		}
		return @result;
	}
	elsif ( $method =~ m/get_distinct_([^_]+)s/ ) {
		my $field = $1;
		my $self = shift;
		my %result;
		for my $row ( $self->as_2d ) {
			my $value = $row->[$keys{$field}];
			$result{ $value } = 1 if $value;
		}
		return sort { $a cmp $b } keys %result;
	}
	elsif ( $method =~ /^[A-Z]+?$/) {
		# do nothing
		return;
	}
	else {
		throw 'UnknownMethod' => $AUTOLOAD;
	}
}
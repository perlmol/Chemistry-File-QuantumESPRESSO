package Chemistry::File::QuantumESPRESSO;

# VERSION
# $Id$

use base 'Chemistry::File';
use Chemistry::Mol;
use Math::VectorReal qw(vector);

use strict;
use warnings;

sub parse_string {
    my ($self, $s, %opts) = @_;

    my $mol_class  = $opts{mol_class}  || 'Chemistry::Mol';
    my $atom_class = $opts{atom_class} || $mol_class->atom_class;
    my $bond_class = $opts{bond_class} || $mol_class->bond_class;
    local $_;

    my $mol = $mol_class->new;

    my @atoms;
    my @cell_vectors;

    my @lines = split "\n", $s;
    while (@lines) {
        my $line = shift @lines;
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;

        if(      $line =~ /^ATOMIC_POSITIONS crystal$/ ) {
            my( $symbol, @coords ) = split /\s+/, $line;
            push @atoms, [ $symbol, vector( @coords ) ];
        } elsif( $line =~ /^CELL_PARAMETERS angstrom$/ ) {
            @cell_vectors = map { vector( split /\s+/, $_ ) }
                                ( shift @lines, shift @lines, shift @lines );
        }
    }

    for my $atom (@atoms) {
        $mol->new_atom( symbol => $atom->[0], coords => $atom->[1] );
    }

    return ( $mol );
}

1;

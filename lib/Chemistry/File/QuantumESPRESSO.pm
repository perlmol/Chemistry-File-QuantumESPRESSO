package Chemistry::File::QuantumESPRESSO;

# VERSION
# $Id$

use base 'Chemistry::File';
use Chemistry::Mol;
use Math::MatrixReal;
use Math::Trig;
use Math::VectorReal qw(vector);

use strict;
use warnings;

Chemistry::Mol->register_format(scf.in => __PACKAGE__);

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

    if( @cell_vectors ) {
        my @lengths = map { $_->length } @cell_vectors;
        my $alpha = acos( ($cell_vectors[1] . $cell_vectors[2]) / ($lengths[1] * $lengths[2]) );
        my $beta  = acos( ($cell_vectors[0] . $cell_vectors[2]) / ($lengths[0] * $lengths[2]) );
        my $gamma = acos( ($cell_vectors[0] . $cell_vectors[1]) / ($lengths[0] * $lengths[1]) );
    }

    for my $atom (@atoms) {
        $mol->new_atom( symbol => $atom->[0], coords => $atom->[1] );
    }

    return ( $mol );
}

sub name_is {
    my ($self, $fname) = @_;
    $fname =~ /\.scf\.in(\.txt)?$/i;
}

sub file_is {
    my ($self, $fname) = @_;
    $fname =~ /\.scf\.in(\.txt)?$/i;
}

# Taken from cod-tools, rev. 9257
sub _symop_ortho_from_fract
{
    my ($a, $b, $c, $alpha, $beta, $gamma) = @_;
    my ($ca, $cb, $cg) = map {cos} ($alpha, $beta, $gamma);
    my $sg = sin($gamma);

    return new_from_rows(
        [ $a, $b*$cg, $c*$cb               ],
        [  0, $b*$sg, $c*($ca-$cb*$cg)/$sg ],
        [  0,      0, $c*sqrt($sg*$sg-$ca*$ca-$cb*$cb+2*$ca*$cb*$cg)/$sg ]
    );
}

1;

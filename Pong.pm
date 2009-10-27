package JCMT::TCS::Pong;

=head1 NAME

JCMT::TCS::Pong - calculate pong pattern properties

=head1 SYNOPSIS

  use JCMT::TCS::Pong;

  $duration = JCMT::TCS::Pong::get_pong_dur( %scan );

=head1 DESCRIPTION

Routines to calculate properties of pong patterns.

=cut

use 5.006;
use strict;
use warnings;
use Carp;
use Scalar::Util qw/ looks_like_number /;

use vars qw/ $VERSION @EXPORT_OK /;
$VERSION = '0.1';

use base qw(Exporter DynaLoader);

@EXPORT_OK = qw( get_pong_dur );

bootstrap JCMT::TCS::Pong;

=head1 FUNCTIONS

=over 4

=item B<get_pong_dur>

  $duration = get_pong_dur( %maparea, %scan );

Where %scan matches the interface used by the JAC::OCS::Config::TCS::obsArea
scan() method and %maparea matches the interface returned by the 
obsArea maparea() method. Since this is a list the order of the hashes
does not matter.

Duration is returned in seconds.

=cut

sub get_pong_dur {
  my %scan = @_;

  croak "Not a PONG pattern. Is a $scan{PATTERN}\n"
    unless $scan{PATTERN} =~ /LISS|PONG/;

  for my $k (qw/ HEIGHT WIDTH VELOCITY DY /) {
    croak "Must provide value for $k"
      unless defined $scan{$k};
    croak "Item $k must be a number. Currently is '$scan{$k}'"
      unless looks_like_number( $scan{$k} );
  }

  my $height = $scan{HEIGHT};
  my $width = $scan{WIDTH};
  my $velocity = $scan{VELOCITY};
  my $dy = $scan{DY};

  my $type = 'CURVY';
  if ($scan{PATTERN} =~ /^(SQUARE|CURVY|ROUNDED)_PONG/) {
    $type =  $1;
  } elsif ($scan{PATTERN} =~ /^(LISS)/) {
    $type = "CURVY";
  } else {
    croak "Unrecognized pong pattern '$scan{PATTERN}'";
  }
  
  # Call the XS routine
  my $duration = _get_pong_dur($height,
                               $width,
                               $dy,
                               $velocity,
                               $type);

  return $duration;
}

=back


=head1 AUTHORS

Tim Jenness E<lt>t.jenness@jach.hawaii.eduE<gt> for the perl
interface.

Russell Kackley wrote the Pong C code.

Copyright (C) 2008 Science and Technology Facilties Council.
All Rights Reserved.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful,but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place,Suite 330, Boston, MA  02111-1307, USA

=cut


1;

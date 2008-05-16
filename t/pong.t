#!perl

use strict;
use warnings;
use Test::More tests => 3;

require_ok( "JCMT::TCS::Pong" );

is(JCMT::TCS::Pong::get_pong_dur( HEIGHT => 320,
                                  WIDTH => 400,
                                  DY => 50,
                                  VELOCITY => 100,
                                  PATTERN => 'CURVY_PONG'), 60,
    "Curvy pong 320x400 at 100 arcsec/sec");


is(JCMT::TCS::Pong::get_pong_dur( HEIGHT => 320,
                                  WIDTH => 400,
                                  DY => 50,
                                  VELOCITY => 600,
                                  PATTERN => 'SQUARE_PONG'), 10,
    "Square pong 320x400 at 600 arcsec/sec");


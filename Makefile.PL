use ExtUtils::MakeMaker;
WriteMakefile(
	      'NAME'      => 'JCMT::TCS::Pong',
	      'VERSION'   => '0.1',
	      'OBJECT'    => 'Pong.o ocsPtcsComputePong.o',
	      'DEFINE'    => '-DSTATUS__OK=0',
              'EXE_FILES' => [qw/ pongtime /],
    );

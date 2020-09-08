#!/usr/bin/env perl

use GraphViz::Makefile;
my $gm = GraphViz::Makefile->new(undef, "Makefile");
$gm->generate("all"); # or another makefile target
open my $ofh, ">", "makefile.png" or die $!;
binmode $ofh;
print $ofh $gm->GraphViz->as_png;
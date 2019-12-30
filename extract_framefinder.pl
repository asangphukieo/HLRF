#!/usr/bin/perl -w

use strict;
use warnings;

my @file_data=();
my $line='';
my $ffname='';
$ffname =  $ARGV[0];

#*** Open input file***
open(FILE,$ffname)or die("Couldn't open Input file\n");
@file_data = <FILE>;
close(FILE);

#print join("\t", "QueryID", "CDSLength", "Score", "Used", "Strict"), "\n";

# >gi|344283|emb|A01270.1| T.ovis mRNA for 45W antigen (partial) \\
# [framefinder (0,702) score=55.32 used=75.65% {forward,local} ]


 foreach $line (@file_data)
     {

if ($line =~ /^>/){

    if (my ($start, $end, $score, $used, $type)
	= ($line =~ /framefinder \((\d+),(\d+)\) score=(\S+) used=(\S+)% \{forward,(\w+)\} /)) {
      
      print join("\t",$score, $used);
   }

}}


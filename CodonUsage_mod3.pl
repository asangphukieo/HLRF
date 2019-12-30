#!/usr/ucb/perl -w

##################################################################
#
# Usage:
# CodonUsage.pl  [options]  FASTA-File(s)
#
# Options:  -c       Check all input sequences are ORFs
#           -s name  Specify score file name
#
# Read FASTA input file(s)
# Count in-frame codon frequencies (without start ATG)
# Count out-of-frame tripletts
# Compute triplett scores for coding/noncoding preference
#
##################################################################
use strict;
use Getopt::Std;

#################################################################

my $prg;
my $check=0;
my $minlength=-1;
my $totalseq=0;
my %ccount;
my %cfreq;
my $ctotal;
my %ncount;
my %nfreq;
my $ntotal;
my %score;
my $totalscore;
my $totalbitscore;
#################################################################


sub info
  {
    my $text=shift;
    print STDERR "$prg: Error: $text\n";
  }


sub usage
  {
    print STDERR <<ENDOFUSAGE
Usage:
CodonUsage.pl  [options]  FASTA-File(s)

Options:  -c       Check all input sequences are ORFs
          -s name  Specify score file name

Read FASTA input file(s)
Count in-frame codon frequencies (without start ATG)
Count out-of-frame tripletts
Compute triplett scores for coding/noncoding preference

ENDOFUSAGE
  }
#####################################################################

# sub process ($header, $seq)
sub process
  {
    my $h=shift;
    my $s=shift;
    my $l;
    my $ww;
    my $start;
    my $codon;

    $s=uc($s);
    ($ww=$s) =~ s/[ACGT]//g;
    if ($check && ($l=length($ww)))
      { info("$l non-ACGT characters found in [$h]\n");  }
    if ($check && !(substr($s,0,3) eq 'ATG'))
      { info("Sequence [$h] does not start with ATG");  }
    $l=length($s);
    if ($minlength<0 || $l<$minlength) {$minlength=$l;}
    if ($check && ($l % 3))
      { info("Sequence [$h]'s length $l not divisible by three!"); }
    for ($start=1; $start<$l-2; $start++)  {
      $codon=substr($s,$start,3);
      if ($start%3==0) {
        if ($check && ($codon =~ /(TGA|TAA|TAG)/))
          { info("Stop Codon $1=$codon in sequence [$h] at $start of $l"); }
        $ctotal++;
        $ccount{$codon}++;
      }
      else { # Out of frame
        $ntotal++;
        $ncount{$codon}++;
      }
    }
    $totalseq++;
  }
#################################################################

# Option parsing
my %options;
$prg = $0;
if ($prg =~ /^(.*)\/(.*)\.pl/) { $prg=$2; }
elsif ($prg =~ /(.*)\.pl/) { $prg=$1; }
getopts("cs:", \%options);
while (($ARGV[0]) && ($ARGV[0] =~ /^-/)) { shift; }
$check = $options{c};
my $scorefile = $options{s};

# Initialize codon counts, freqs to zero
%ccount = ();
%cfreq  = ();
$ctotal=0;
%ncount = ();
%cfreq  = ();
$ntotal=0;
$totalscore=0;

# Processing loop
my $header = "NIX";
my $line;
my $seq = "";
while (! ($header =~ /^>/) )  { $header = <>; }
chomp($header);
while(<>)
  {
    chomp; $line=$_;
    if ($line =~ /^>/)
      { # Process $seq
        process($header,$seq);
        $header = $line; $seq="";
      }
    else   {  $seq .= $line;   }
  }
process($header, $seq);


# Write Codon Usage Statistics
my $codon;
#printf("Cod: InFram  freq   || OOFram  freq   || ScoreBits\n");
foreach $codon (sort keys %ncount) {
  $nfreq{$codon}=$ncount{$codon}/$ntotal;
  if (exists ($ccount{$codon})) {
    $cfreq{$codon}=$ccount{$codon}/$ctotal;
    $score{$codon}=log($cfreq{$codon}/$nfreq{$codon})/log(2);
  }
  else {
    $ccount{$codon}=0;
    $cfreq{$codon}=0;
    $score{$codon}=-10;
  }
  
  $totalscore = $totalscore + $score{$codon};
}
#printf("---: %6d  ------ || %6d  ------ || ---------\n",$ctotal, $ntotal);
$totalbitscore=log($ctotal/$ntotal)/log(2);
printf("%12.8f ",$totalbitscore);
 
printf("%9.5f ",$totalscore);

# Write additional scorefile if requested
if ($scorefile) {
  open SF,">$scorefile";
  foreach $codon (sort keys %ncount) {
    printf SF "$codon %12.8f\n", $score{$codon};
  }
  close SF;
}


# Show statistics
#print STDERR "Done. Processed $totalseq input sequences. Minlength=$minlength\n";
#################################################################

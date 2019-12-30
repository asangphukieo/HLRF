#!/usr/bin/perl 

use lib "/home/marasri/NongDoc/RNA_jeab/HeteroMirPred/HeteroMirPred/ViennaRNA-1.8.5/Perl/blib/arch";
use lib "/home/marasri/NongDoc/RNA_jeab/HeteroMirPred/HeteroMirPred/ViennaRNA-1.8.5/Perl/blib/lib";

no warnings;
use strict;
use Getopt::Long;
use RNA;

############################################################################
# Global Parameters and initialization if any.
############################################################################

my $inFile="&STDIN";
my $outFile="&STDOUT";

#Define the monomers and dimers
my %gl_monomers = ('A' => 0,'C' => 0, 'G' => 0, 'U' => 0);
my %gl_dimers = ('AA' => 0, 'AC' => 0, 'AG' => 0, 'AU' => 0,
                 'CA' => 0, 'CC' => 0, 'CG' => 0, 'CU' => 0,
                 'GA' => 0, 'GC' => 0, 'GG' => 0, 'GU' => 0,
                 'UA' => 0, 'UC' => 0, 'UG' => 0, 'UU' => 0);
my %gl_4mers = ('GUUC' => 0, 'GGUU' => 0, 'GCAU' => 0, 'GAUA' => 0, 'CUAC' => 0, 'CAGU' => 0, 'AGGA' => 0,
                 'AAAA' => 0, 'AUGA' => 0, 'CUGA' => 0, 'UACA' => 0, 'CAAC' => 0, 'CGGA' => 0, 'GAAG' => 0,
                 'UCGU' => 0, 'CGUU' => 0, 'UGAU' => 0, 'UCCG' => 0, 'UAAG' => 0 ,  'UUUU' => 0 );
my $numseqs = 0;
############################################################################
# File IO
# Parse the command line.
############################################################################
Getopt::Long::Configure ('bundling');
GetOptions (
	'i|input_file=s' => \$inFile, 
	'o|output_file=s' => \$outFile
);

if(scalar(@ARGV) == 1 || !defined($inFile) || !defined($outFile)) { 
	die ("USAGE: $0 -i <input file> -o <output file>\n");
}

open (INFILE, "<$inFile") or die( "Cannot open input file $inFile: $!" );
open (OUTFILE, ">$outFile") or die ("Cannot open output file $outFile: $!");


#print (OUTFILE map { "$_\t" } (sort keys(%gl_4mers)));

# Read line by line.
while (my $line = uc(<INFILE>)) {

	chomp($line);
	$line =~ s/T/U/g;
	
	# Fasta First Line
    if ($line =~ m/^>/) { }
    
    # Fasta Second Line i.e. RNA sequence
    elsif ($line =~ m/^[AaCcUuGg]/) {	    

		#Absolute Values
		my %aw_4mers = %gl_4mers;		

		$numseqs++;

		#remove white space etc
		$line =~ s/[^AaCcUuGg]//g;

		my $seqLen = length($line);
		#print(OUTFILE "\n");
		
		#compute monomer and dimer distribution
		for my $i (0..$seqLen-1) {		
			my $mer = substr($line, $i, 4);
			$aw_4mers{$mer}++ if defined $aw_4mers{$mer};			
		}	
		
		#Print Absolute Values	
		#foreach my $mer (sort (keys(%aw_4mers))){
			#print(OUTFILE "$aw_4mers{$mer}\t");	
		#}
		
                foreach my $mer (sort (keys(%aw_4mers))){
			printf(OUTFILE "%.2f\t", $aw_4mers{$mer}/($seqLen-1)*100);
		}
		
    }
	
    else { }
  
}#end of while loop

print (OUTFILE "\n");
close (INFILE) or die( "Cannot close input file $inFile: $!" );
close (OUTFILE) or die( "Cannot close output file $outFile: $!");
exit;



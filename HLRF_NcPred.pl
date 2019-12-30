

use strict;
use warnings;

#*** Initialization ***
   my @file_data=();
   my $filename='';
   my @file=();
   my @sentences=();
    my $line='';
    my $tempfile='';
    my $temp='';
    my $tempdata='';
    my $tempdata2='';
    my $empty='';   
    my $sentence='';
    my $fastaSeq='';
   my @readspec=();
   my $foldfile='';
   my $allstatfile='';
   my $num = 0 ;

my $statfilename='';
my $feature20name='';
my $SCfilename='';
my $spectralfile='';
my $mfe14file = '';
my $mfe23file ='';
my $BP1file ='';
my $BP2file ='';
my $line1='';
my $line2='';
my $line3='';
my $line4='';
my $line5='';
my $line6='';
my $line7='';
my $line8='';
my $line9='';
my $spec='';
my $tab = "\t";
my $count_seq= 1 ;

##**Change multiple line sequence to one line fasta format which is the format for features extraction module***
#system("perl fasta1line.pl $ARGV[0] $ARGV[1]");

#system("rm $ARGV[0]");

#Edit_byNong
system("rm *.ps");
system("rm *.fold");

$filename =  $ARGV[0];
#*** Open input file***
open(FILE,$filename)or die("Couldn't open Input file\n");
@file_data = <FILE>;
close(FILE);
system("rm $filename.csv");system("rm $filename.SeqList.csv");
#system("rm $filename.20features.csv");

open(FILE3, ">> $filename.csv");
$tempdata = "mfe,Prob,efe,mfe1,mfe2,dG,dP,dQ,dF,zG,zP,zQ,mfe3,mfe4,nefe,freq,div,diff,dH,dHL,dH/Loop,dS,dSL,dS/loop,Tm,TmL,Tm/Loop,au/L,gc/L,gu/L,P_au/bp,P_gc/bp,P_gu/bp,tot_bp/Loops,mfe5,SCI,SC/non_tot,SCxMFE/Mean,SCxabsZG,SCxdP,SC/len,SC/auLgcLguL,SC/nonauLgcLguL,SC/1dP,SC/nonA,SC/nonU,score1,score2,score3,score4,score5,score6,score7,score8,score9,score10,score11,score12,score13,score14,score15,score16,score17,score18,score19,score20,score21,score22,score23,score24,score25,score26,score27,score28,score29,score30,score31,score32,score33,score34,score35,score36,score37,score38,score39,score40,score41,score42,score43,score44,score45,score46,score47,score48,score49,score50,score51,score52,score53,score54,score55,score56,score57,score58,score59,score60,score61,score62,score63,score64,score65,score66,score67,score68,score69,score70,score71,score72,score73,score74,score75,score76,score77,score78,score79,score80,score81,score82,score83,score84,score85,score86,score87,score88,score89,score90,score91,score92,score93,score94,score95,score96,score97,score98,score99,score100,a1,a2,a3,a4,a5,a6,a7,a8,g1,g2,g3,g4,g5,g6,g7,g8,c1,c2,c3,c4,c5,c6,c7,c8,t1,t2,t3,t4,t5,t6,t7,t8,mfe_a,MeanBP,A,C,G,U,AA,AC,AG,AU,CA,CC,CG,CU,GA,GC,GG,GU,UA,UC,UG,UU,pairprob1,pairprob2,pairprob3,pairprob4,pairprob5,pairprob6,pairprob7,pairprob8,pairprob9,pairprob10,Non_BPP,nonBP_a,nonBP_c,nonBP_g,nonBP_u,logbits,sumbits,Bits,Bits2,CM,New4_2,LR2,AAAA,AGGA,AUGA,CAAC,CAGU,CGGA,CGUU,CUAC,CUGA,GAAG,GAUA,GCAU,GGUU,GUUC,UAAG,UACA,UCCG,UCGU,UGAU,UUUU,Class\n";
print FILE3 $tempdata;

open(FILE4, ">> $filename.SeqList.csv");
$tempdata2 = "SeqName\n";
print FILE4 $tempdata2;

 # print  @file_data;
foreach $line (@file_data)
     {
        if ($line =~ /^>/){
          $temp  .= $line;
          $num++;
	  $fastaSeq=">$filename.fasta";
	  #$foldfile=">>$num.fold";	                   
          }
       
       elsif ($line =~ /^\w/){
             $temp .= $line;
          
	     open(OUTPUT, $fastaSeq);	
             print OUTPUT $temp;
	     
   	     print("\n-calculating features-------------\n");

	     system("RNAfold < $filename.fasta > $filename.fold"); 
             system("perl genRNAStats_mod.pl < $filename.fasta > $filename.stat2");
	     system("perl genRNAStats.pl < $filename.fasta > $filename.data1");
	     system("./RNAtopological.exe < $filename.fold > $filename.spec");
	     
	     system("perl genRandomRNA.pl -n 50 -m d $filename.fasta > $filename.random.fasta");
		
		system("RNAfold < $filename.random.fasta > $filename.random.fold"); 
		#Edited by Nong
		system("rm *_ss.ps");
	     
 	     system("perl genRNARandomStats_mod.pl -n 50 -i $filename.random.fold -o $filename.zdata -m $filename.fold");
	     $tempfile=">temp.txt";
		system("rm *.random.fold");
   		system("rm *.random.fasta");
   		system("java mfe14 $filename.data1");
   		system("java mfe23 $filename.spec");
		system("perl melt2.pl $filename.fasta > $filename.mfold2");
		system("perl melt.pl $filename.fasta > $filename.mfold");
		system("java Mfoldfilter $filename");
		system("python selfcontain.py -i $filename.fasta > $filename.sc");
   		system("java bpcount1_mod $filename");
                #Count Percent of Trimer occurs in LongStem------------------------------- 
                #system("perl trimerstemcount_2.pl < $filename.longstem > $filename.trimer");
   		system("RNAfold -p2 < $filename.fasta > $filename.RNAfold1");  
   		system("java RNAfoldfilter $filename");
		
		#Edit_byNong
   		#system("rm *.ps");
		system("mv *_ss.ps ss_$count_seq.ps");
		system("rm *dp.ps");
		system("rm *dp2.ps");
		
		system("rm $filename.fasta.37.ext");
		system("rm $filename.fasta.ct");
		system("rm $filename.fasta.run");
		system("./randfold -d $filename.fasta 100 > $filename.prob");
                system ("./Conv.exe $filename.fasta $filename.f2 1");
		system ("./Extract_PP.exe -f $filename.f2 -out $filename.Probpair");
                system ("perl CodonUsage_mod3.pl $filename.fasta > $filename.codon");
                system ("perl motifStats.pl < $filename.fasta > $filename.motif");
		system ("./blastn -query $filename.fasta -db DB/Rfam90remove -out $filename.blast -word_size 4");
		system ("./blastn -query $filename.fasta -db DB/database2 -out $filename.blast2 -word_size 5");
                system ("cmscan CM/Rfam_mod.cm $filename.fasta > $filename.OUT");
                system ("java myCMfilter $filename");
		system ("./framefinder -r false -w framefinder.model $filename.fasta > $filename.frame");
                system ("perl extract_framefinder.pl $filename.frame > $filename.ff");
		system("java filter_F2 $filename"); 
                #system("java filter_F4 $filename"); 
		print(" Features for sequence no. $num is Done ----------------------------\n");	
		
		my $statfilename= "$filename.feature";
		
		open(FILE1, $statfilename);
		
                my $feature20name= "$filename.Name";
 		open(FILE2, $feature20name);
			

		while(defined(my $line1 = <FILE1>)) {
			print FILE3 $line1;
			$line2 = <FILE2>;	
			print FILE4 $line2;
			}
              
	         system("rm $filename.fasta");
                 system("rm $filename.feature");
		 system("rm $filename.bp1");
 		 system("rm $filename.data1");
		 system("rm $filename.longstem");
		 system("rm $filename.mfold2");
 		 system("rm $filename.mfold");
                  system("rm $filename.sc");
		 system("rm $filename.RNAfold1");
		 system("rm $filename.RNAfold2");
		 system("rm $filename.spec");
		 system("rm $filename.stat2");
		 system("rm $filename.fasta.37.plot");
		 system("rm $filename.spec.mfe23");
		 system("rm $filename.data1.mfe14");	
		system("rm $filename.fasta.dG");
		
		
        #system("rm $filename.fold");	
		#change file bracket name by Nong
		system("mv $filename.fold bracket_$count_seq.fold");
		
                system("rm $filename.zdata");
		system("rm $filename.prob");
                system("rm $filename.f2");
                system("rm $filename.Probpair");
		system("rm $filename.codon");
		system("rm $filename.blast");
		system("rm $filename.blast2");
		system("rm $filename.motif");
		system("rm $filename.cmout2");
                system("rm $filename.OUT");
		system("rm $filename.frame");
                system("rm $filename.ff");
		system("rm $filename.Name");
		#my $feature20name="$num.predict";
		
           #  print TEMPOUT $num;
             close(OUTPUT);
            # close(TEMPOUT);
             $temp=$empty;
		$count_seq++;	 

           }
	

}
$tempdata="-19.6,0.584158,-22.2,-0.0040825,-0.0233333333,-0.1633,0.2583,0.2315,0.166717,0.0311,0.7233,0.1513,-0.0204166667,-0.6322580645,-0.185,0.0240506,8.515,0.0216666667,-250.7,-2.0891666667,-31.3375,-760,-6.3333333333,-95,56.7,0.4725,7.0875,0.125,0.1,0.0333333333,48.3870967742,38.7096774194,12.9032258065,3.875,-0.0158064516,0.4140336134,0.0046520631,0.4117563487,0.0128764454,0.1069448824,0.0034502801,1.6027107617,0.5582475687,0.8565031308,0.6533270638,0.8889228878,35.53696,36.3264,36.657152,47.716878,34.674204,40.318166,34.59723,39.378919,42.811451,39.068928,36.834816,37.895986,35.023304,34.764309,36.023305,42.203781,47.232692,38.829476,38.46575,46.776136,33.458824,35.357043,39.631075,46.73865,27.097005,37.686742,38.469862,46.267872,38.830565,41.044592,42.856231,38.173877,40.702443,40.955017,39.006289,44.368264,39.48582,38.107944,44.418494,39.999069,41.596226,36.031604,35.203103,43.853799,41.270768,30.051437,40.842589,48.861899,36.881724,37.056789,32.71773,43.403871,38.327873,39.156765,54.786654,41.713615,36.549352,44.045128,55.116909,37.599442,48.619825,36.836579,37.340683,46.12902,42.805008,51.285542,40.841757,45.080719,53.302859,47.106681,41.198915,41.669963,45.063066,42.165801,42.348565,48.962075,41.260655,39.83551,40.30035,41.8786,48.274129,42.266139,41.301208,56.179555,48.636036,54.961611,42.118665,46.22709,44.071452,39.020731,38.970559,39.83883,44.080076,41.691729,39.745321,47.694275,39.491091,41.913148,44.209385,42.440258,0.248994,0.160608,0.044428,0.179701,0.070862,0.005603,0.087502,0.202302,0.283865,0.158989,0.060149,0.046033,0.114441,0.073037,0.109969,0.153515,0.144169,0.013227,0.061371,0.015288,0.206656,0.033249,0.129474,0.396565,0.217603,0.099026,0.020536,0.128606,0.128667,0.034015,0.113786,0.257763,-26.42,14.098451,0.325,0.225,0.175,0.275,0.117647,0.084034,0.033613,0.092437,0.07563,0.02521,0.067227,0.058824,0.042017,0.05042,0.02521,0.05042,0.092437,0.067227,0.042017,0.07563,7.919623,3.613291,2.48429,6.025953,2.228254,2.168801,1.970621,1.971774,5.430009,3.755826,0.50324,0.633731,0.549037,0.260901,0.46577,-1,-268.16993,30.7,75,6.9,151.2386100672,76.792660328,1.68,0,0.84,0.84,0,0.84,0,0.84,0.84,0,0,0,0,0.84,0,0,0,1.68,0,1.68,Other";
print FILE3 $tempdata;

system("java -classpath \$CLASSPATH:weka.jar:libsvm.jar:mysql-connector-java-5.0.8-bin.jar weka.classifiers.meta.Vote -l HRF.model -T $filename.csv -p 36 > $filename.out");

system("clear");

print("\n\nInput file: $filename; There are $num query sequences in the input file. \n");
print("\n.... Completed Extracting All Features ....\n");
print("\nAll features extracted are written to the file: $filename.csv.");

print("\n\n########## Prediction results by Hybrid Random Forest Ensemble (HLRF) ##########\n\n");


system("python read_NCoutput_to_file.py $filename");

exit;

#!/usr/bin/perl
use strict;
use warnings;
#
#Usage on multiple samples: for i in $(ls s_*_1_sequence.txt); do Adaptertrim.pl -i $i -Bm -U - -o /home/terpstramm/Documents/RNA-seq/trim/${i}_adaptertrimmed.txt -a AGATCGGAAGAGCAC|fastx_trimmer -l 86  -m 25 -v> /home/terpstramm/Documents/RNA-seq/trim/${i}_14nt_trimmed.txt; cat /home/terpstramm/Documents/RNA-seq/trim/${i}_adaptertrimmed.txt /home/terpstramm/Documents/RNA-seq/trim/${i}_14nt_trimmed.txt > /home/terpstramm/Documents/RNA-seq/trim/${i}_final_trimmed.txt; done; 
#
#Warning: using this script multiple times on the same pipe scrambles the statistics output if redirected to the same file
my $use = <<"END1";
Use $0 -i FILE -a ADAPTER -[BA]
	trims on adapter in a simple way using a case insensitive index match. This is used for small RNA sequencing with longer sequencing reads (>70nt).
required switches
	-i FILE	input file ('-i -' for STDIN)
	-a	adapter sequence
	-[BA]	trim direction: returns everything Before/After the adapter
		(Before is the usual suspect) 
optional switches
	-f FORMAT specify input format. valid options for FORMAT are 'fasta' and 'fastq'. 'fastq' is default
	-o specify output file ('-o -' for STDOUT) undef = outputfile generated on 
		 input name and trimming options
	-u include untrimmed sequences to output file
	-U specify separate output file for untrimmed sequences
		 overrides -u writes to file ('-U -' for STDOUT) dies if -p 
		 specified
	-m Allows one mismatch in the adapter by adapter permutation and matching all the adapters. this slows down the script considerably.
	-N doesn't allow N containing sequences to be returned
	-v verbose results (to STDERR)
	-p keep PE structure intact, adapter only is outputted as 
		 seq:N qual:B sequences
	-l INT	largest read size to allow in reads. if larger a 
		 substring is taken [default 1000]
	-s INT	smallest read size to allow in reads. if smaller read is ignored
		 or returned as seq:N qual:B sequences [default 1]

examples:
	$0 -i FILE -a tacga -A
		FILE:
			\@head
			TTTTACGANNN
			+head
			BBBBBBBBBBB
	returns 
		FILE_tac:
			\@head
			NNN
			+head
			BBB
	$0 -B -a tac -i FILE  
		FILE:
			\@head
			TTTTACGANNN
			+head
			BBBBBBBBBBB
	returns 
		FILE_tac:
			\@head
			TTT
			+head
			BBB
END1
#######################PSEUDO
# permutate adapter
# use the permutated results to add the different nucleotides to the permutation
#open file
#foreach entry in file{
# match adapter permutations
# if adapter is found this will be trimmed off
# optional filtering for N containing sequences
# other timming steps
# write output/untrimmed to file/STDOUT
#}
#print optional rapport 2 STDERR
#
#no quality check/advanced length filtering/simple trimming 
#
##########################cmdOptparse

if (scalar(@ARGV)==0){
	print $use;
	exit(0);
}

use Getopt::Std;
my %cmdOptions;
getopts('i:a:o:U:f:BAmNuvps:l:', \%cmdOptions);
#defaults
#($cmdOptions{'s'}= 1)if(not($cmdOptions{'s'}));
#($cmdOptions{'l'}= 1000)if(not($cmdOptions{'l'}));#should be inf
#sanity#########################################################################
die "invalid size cuttoffs lower size limit is higher then upper size limit\n"if($cmdOptions{'s'} && $cmdOptions{'l'} &&$cmdOptions{'s'}>$cmdOptions{'l'});
die "invalid lower size limit <0\n"if($cmdOptions{'s'} && $cmdOptions{'s'}<0);
die "invalid upper size limit <0\n"if($cmdOptions{'l'} && $cmdOptions{'l'}<0);


################################################################################

my $adapter	= $cmdOptions{'a'};
my $adapter_len	= length($adapter);

my $mismatches;
if ($cmdOptions{'m'}){
$mismatches	= 1;#cant handle more than 1 not implemented
}else{
$mismatches	= 0;
}
my $format;
if ($cmdOptions{'f'} &&($cmdOptions{'f'} eq 'fasta' || $cmdOptions{'f'} eq 'fastq')){
	$format	= $cmdOptions{'f'};
}elsif($cmdOptions{'f'}){
	die "invalid format specified with -f '$cmdOptions{'f'}' must be 'fasta' or 'fastq'\n$!\n";
}else{
$format	= 'fastq';#default;
}
my $in;
if($cmdOptions{'i'}){
	if($cmdOptions{'i'} eq '-'){
		$in = *STDIN;
	}else{
		open($in,'<',$cmdOptions{'i'}) or die "$0 open error: $!";
	}
}else{
print STDERR "no input file \n";
exit(0);
}

my $out;

if(($cmdOptions{'o'})&&($cmdOptions{'o'} eq '-')){
	$out = *STDOUT;
}elsif($cmdOptions{'o'}){
	open($out,'>',$cmdOptions{'o'}) or die "$0 write error $!";
}else{
	if($cmdOptions{'B'}){
		open($out,'>',$cmdOptions{'i'}."_B_$adapter")or die "$0 write error $!";
	}elsif($cmdOptions{'A'}){
		open($out,'>',$cmdOptions{'i'}."_A_$adapter")or die "$0 write error $!";
	}else {
		print STDERR "$0 error: no trim direction specified (specify by using -A or -B])";
		exit(0);
	}
}
my $outU;
if((defined$cmdOptions{'U'})&&($cmdOptions{'U'} eq '-')){
	$outU = *STDOUT;
	if($outU eq $out){
		print STDERR '-U equals -o. Specify inf only -o and -u without outfile is faster.'."\n";
	}
	$cmdOptions{'u'} = 1;
}elsif($cmdOptions{'U'}){
	open($outU,'>',$cmdOptions{'U'}) or die "$0 write -U(untrimmed) error $!";
	$cmdOptions{'u'} = 1;
}
##############################main
#firstLine error with reading fasta....without this
if($format eq 'fasta'){
	{
		#######################tekst block use with ReadsFasta( &ReadsFasta() and sub ReadsFasta)
		#######################or get annoying errors
		local $/ = ">";
		$_ = <$in>; #let's ignore that annoying first line;
		$_ = '';#gone
	}
}

(my @adapterpermut = &permutAda($adapter,$mismatches))if($cmdOptions{'a'});#do permutations
(my @Ada_array = &generateAll('ATCGN',(@adapterpermut)))if($cmdOptions{'a'});#generate all possible strings by repacing the N's to [ATGCN}
my %seqstats;

while (not eof($in)){
	my $fqHeader1;
	my $fqseq; #set to $fqseq ="N" instead of empty string   if not passing QC for observing the effect of the trimming strategy
	my $fqqual;#set to $fqqual="B" instead of empty string 
	#read in data
	(($fqHeader1, $fqseq, undef, $fqqual) = &read_fastq(*$in))if($format eq 'fastq');
	($fqHeader1, $fqseq) = &read_fasta(*$in)if($format eq 'fasta');
	#			print STDERR join('	',($fqHeader1, $fqseq))."\n";
	$seqstats{'ReadsInInput'}++;
	
	my $bool = 0;
	my $i = 0;
	
	my $len_Ada_array = scalar(@Ada_array);
	my $index = -1;
	if($cmdOptions{'a'}){
		while($bool == 0){#permuted adapter match
			$index = index(uc($fqseq),uc($Ada_array[$i]));
			$i++;
			$bool++ if($i == $len_Ada_array || $index > -1);
		}
	}
	#count if untrimmed because no adapter found
	if($index < 0){
		$seqstats{'ReadsWithAdapterNotFound'}++;
	}

	
	if(((length($fqseq) <= ($index+$adapter_len)) && $cmdOptions{'A'})||((($index) == 0) && $cmdOptions{'B'})){##Error catching for adapter only reads
		###########my mistake ignore plz
		#if($cmdOptions{'p'}){
		#	$fqseq = 'N';
		#	if(defined($fqqual)){
		#		$fqqual = 'B';
		#	}
		#}else{
			$fqseq = 'N';
			if(defined($fqqual)){
				$fqqual = 'B';
			}
		#}
		$seqstats{'ReadsOmittedOnlyAdapter'}++;
		
	}elsif($cmdOptions{'A'} && $index >= 0){#start trim should be done when the adapter index is lt 0
		#$seqstats{'ReadsWithAdapterTrimmed'}++;
		$fqseq = substr($fqseq,($index+$adapter_len));
		if(defined($fqqual)){
			$fqqual = substr($fqqual,($index+$adapter_len));
		}
			
	}elsif($cmdOptions{'B'} && $index > 0){#end trim 
		#$seqstats{'ReadsWithAdapterTrimmed'}++;
		$fqseq = substr($fqseq,0,$index);
		if(defined($fqqual)){
			$fqqual = substr($fqqual,0,$index);
		}
	}
	
	#select for min/max size
	
	my $tmp = 0;
	
	(($fqHeader1, $fqseq, $fqqual,$tmp) = &selectMaxLen($fqHeader1, $fqseq, $fqqual, $cmdOptions{'l'}))if($cmdOptions{'l'} && $fqseq ne 'N');
	
	$seqstats{'ReadsWithTrimmedMaxLen'}+=$tmp;
	
	
	(($fqHeader1, $fqseq, $fqqual,$tmp) = &selectMinLen($fqHeader1, $fqseq, $fqqual, $cmdOptions{'s'}))if($cmdOptions{'s'} && $fqseq ne 'N');
	
	$seqstats{'ReadsOmittedMinLen'}+=$tmp;
	
	#remove N containing 
	if($cmdOptions{'N'}){
		if($fqseq =~ /[nN]/ && $fqseq ne 'N'){
			$seqstats{'ReadsOmittedNContaining'}++;
			$fqseq = 'N';
			$fqqual = 'B';
		}
	}
	

	#print output
	if($cmdOptions{'p'}){
		
		if($cmdOptions{'u'}){
			if($cmdOptions{'U'} && $index < 0){
				die "cannot logically keep PE structure intact in -U FILE! remove -p or specify -u istead of -U FILE\n$!";
			}else{
				$seqstats{'ReadsInOutput'}++;
				if(defined($fqqual)){
					print {$out} "\@$fqHeader1\n$fqseq\n+$fqHeader1\n$fqqual\n";
				}else{
					print {$out} "\>$fqHeader1\n$fqseq\n";
				}
			}	
		}else{
			$seqstats{'ReadsInOutput'}++;
			if(defined($fqqual) ){
					print {$out} "\@$fqHeader1\n$fqseq\n+$fqHeader1\n$fqqual\n";
			}else{
				print {$out} "\>$fqHeader1\n$fqseq\n";
			}
		}
		
	}elsif($fqseq ne 'N'){
		if($cmdOptions{'u'}){
			if($cmdOptions{'U'} && $index < 0){
				$seqstats{'ReadsInOutput'}++;
				if(defined($fqqual)){
					print {$outU} "\@$fqHeader1\n$fqseq\n+$fqHeader1\n$fqqual\n";
				}else{
					print {$outU} "\>$fqHeader1\n$fqseq\n";
				}
			}elsif($cmdOptions{'U'} && $index >= 0){
				$seqstats{'ReadsInOutput'}++;
				if(defined($fqqual)){
					print {$out} "\@$fqHeader1\n$fqseq\n+$fqHeader1\n$fqqual\n";
				}else{
					print {$out} "\>$fqHeader1\n$fqseq\n";
				}
			}else{
				$seqstats{'ReadsInOutput'}++;
				if(defined($fqqual)){
					print {$out} "\@$fqHeader1\n$fqseq\n+$fqHeader1\n$fqqual\n";
				}else{
					print {$out} "\>$fqHeader1\n$fqseq\n";
				}
			}	
		}elsif($index > 0){
			$seqstats{'ReadsInOutput'}++;
			if(defined($fqqual)){
				print {$out} "\@$fqHeader1\n$fqseq\n+$fqHeader1\n$fqqual\n";
			}else{
				print {$out} "\>$fqHeader1\n$fqseq\n";
			}			
		}
	}
}
close($in);
close($out);
if($cmdOptions{'v'}){
	print STDERR '#'.$cmdOptions{'i'}."\n";
foreach my $key (sort(keys(%seqstats))){
	print STDERR $key."\t".$seqstats{$key}."\n";
}
print STDERR "\n";
}
exit(0);
#######################subs

sub read_fastq {
	my $fileh = shift @_;
	my $read_fqHeader1 = <$fileh>;
	chomp($read_fqHeader1);
	$read_fqHeader1 = substr($read_fqHeader1, 1);
	my $read_fqseq = <$fileh>;
	chomp($read_fqseq);
	my $read_fqHeader2 = <$fileh>;
	chomp($read_fqHeader2);
	$read_fqHeader2 = substr($read_fqHeader2, 1);
	my $read_fqqual = <$fileh>;
	chomp($read_fqqual);
	if ($read_fqHeader1 ne $read_fqHeader2 && $read_fqHeader2 ne ''){
		warn "$0 header error: '$read_fqHeader1' doesn't match '$read_fqHeader2'\n";
		exit(0);
	}
	return($read_fqHeader1, $read_fqseq, $read_fqHeader2, $read_fqqual);
}
#not optimized
sub read_fasta {
	my $filehandle = $_[0];	
	local $/ = '>';
	my $line;
	$line = <$filehandle>;

	chomp($line);
	$line =~ s/^(\S+)//;
	my $faHeader = $1;	
	$line =~ s/\n|\r//gm;#Windows.....
	my $sequence = $line;
	my @fasta = ($faHeader,$sequence);
	return @fasta;
}


sub permutAda {
#input (Adapterseq,amountofmismatches[0..1])
	#amountofmismatches[0..1] -> larger amounts gave me a headace so they aren't included (and besides if the adapter contained that many mm what does that say of the read before the adapter).
	my $permutAda_Adapter = shift @_;
	my $permutAda_length = length($permutAda_Adapter);
	my $permutAda_mismatches = shift @_;
	#my $permutAda_gaps = shift @_;
	my @permutAda_adaPermutations;
	push(@permutAda_adaPermutations,$permutAda_Adapter);
	#permut mismatches
	my $i_mismatches;
	for($i_mismatches=0; $i_mismatches < $permutAda_mismatches; $i_mismatches++){
		my $i_ada;
		for($i_ada=1; $i_ada <= $permutAda_length; $i_ada++){
			my $begin = substr($permutAda_Adapter,0,$i_ada);
			chop($begin);
			$begin = ${begin}."N";
			if(($i_ada - $permutAda_length)){
				my $end = substr($permutAda_Adapter,($i_ada - $permutAda_length));
				push(@permutAda_adaPermutations,$begin . $end);
			}else{
				push(@permutAda_adaPermutations,$begin);
			}
		}
	}
	return @permutAda_adaPermutations;
}

sub generateAll{
#subtitutes N and - for indicated bases
#input bases, ADAPTER and sequences liek: ("ATCGN",$ADAPTER,@Array)
	my $rep = shift @_;
	my $rep_length = length($rep);
	my $adap = shift @_;
	my @Sequences;
	push (@Sequences,$adap);
	foreach my $seq(@_){
		if ($seq =~ /[N-]/){
			my $i_rep;
			for($i_rep=0; $i_rep < $rep_length; $i_rep++){
				my $char = substr($rep,$i_rep,1);
				my $seq_TMP = $seq;
				$seq_TMP =~ s/[N-]/$char/;
				if($seq_TMP =~ /$adap/){
				}else{				
					push(@Sequences,$seq_TMP);
				}
			}
		}else{
		push(@Sequences,$seq);
		}
	}
	return @Sequences;
}
##################################incomplete_subs
sub trimQual{
#
#to be made
#
	my ($fqHeader1, $fqseq, $fqqual, $min_qual) = @_;
	#reverse($fqqual);
	my $i;
	for($i = (length($fqqual)-1);$i<0;$i--){
	#print STDERR (ord(substr($qual,$i,1))-64)."\t";
		last if($min_qual < (ord(substr($fqqual,$i,1))-64));
	}
	$fqseq = substr($fqseq,1,$i);
	$fqqual = substr($fqqual,1,$i);
}
#use ($fqHeader1, $fqseq, $fqqual) = &selLen($fqHeader1, $fqseq, $fqqual,$PE, $minLen, $maxLen);
sub selectMinLen{
#
#Selects for min/max length....to be made/tested
#	
	my $fqHeader1;
	my $fqseq;
	my $fqqual;
	my $minLen;
	($fqHeader1, $fqseq, $fqqual, $minLen) = @_;
	if($fqqual){
		return($fqHeader1,'N','B',1)if(length($fqseq) < $minLen);
		return($fqHeader1,$fqseq,$fqqual,0);
	}else{
		return($fqHeader1,'N',undef,1)if(length($fqseq) < $minLen);
		return($fqHeader1,$fqseq,undef,0);
	}
}
sub selectMaxLen{
#
#Selects for min/max length....to be made/tested
#	
	my $fqHeader1;
	my $fqseq;
	my $fqqual;
	my $maxLen;
	($fqHeader1, $fqseq, $fqqual, $maxLen) = @_;
	if($fqqual){
		return($fqHeader1,substr($fqseq,0,$maxLen),substr($fqqual,0,$maxLen),1)if(length($fqseq) > $maxLen);
		return($fqHeader1,$fqseq,$fqqual,0);
	}else{
		return($fqHeader1,substr($fqseq,0,$maxLen),undef,1)if(length($fqseq) > $maxLen);
		return($fqHeader1,$fqseq,undef,0);
	} 
}

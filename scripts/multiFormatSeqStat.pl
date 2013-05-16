#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Std;
my %Options;
getopts('f:Qr',\%Options);
#default options
$Options{q} = 1;#percentile qual print
$Options{Q} = 0;#no full qual print
$Options{f} = "";#default no force = determine based on reading the first few lines (works for FASTA FASTQ LIST)

#reverse mode for checking correct adapter trimming.
my $reverse = 0;
if($Options{r}){
  $reverse = $Options{r};
}
my $use = <<"END1";
$0 [options] <file|-> 
gets length/base/qual distribution statistics from multi fasta/fastq files. with a partial support for collapsed(like those used in miRDeep2 and miRAnalyzer) formats.

-q	percentile qual print [default:True]
-Q	full qual print [default:False]
-f	Force format to red the file with at the moment FASTA MIRDFASTA FASTQ LIST FTF\\d{6} are supported.[default:disabled] = determine based on reading the first few lines closing the handle then reading the file again (works for FASTA FASTQ LIST). 

note:with force on the correct format the file can be streamed.

outputs to the input file name with .log added
END1
#for future use: program description
my $debug = 0;

my $tofile = 0;
my $out = *STDOUT; # could be anything...
#may the force be with you...:)
my $force = '';#for forcing a format to read the file with.......FASTA FASTQ MIRDFASTA FTF010203
if($Options{f} && $Options{f} =~ /^FASTA$|^MIRDFASTA$|^FASTQ$|^LIST$|^FTF\d{6}$/){
$force = $Options{f};
}
foreach my $file(@ARGV){
	my $format;
	($format = &determineFileFormat($file))if(not($force));
	if($tofile){
		my $fileout = $file;
		$fileout = $fileout.'.log';
		open($out,'>',$fileout);
	}
	#sleep(1);
	($format = $force)if($force);
	print STDERR "$format\n";
	my %data;
	my %nucleotides;
	my $in;
	if($file eq '-'){
	$in = *STDIN;
	print STDERR "Reading from STDIN\n" if($force);
	die "A file type must be declared in the \$force parameter (e.g. FASTA FASTQ MIRDFASTA FTF010203). because STDIN is not suitable for the type of format tests in this programme!\n" if(not($force));
	}else{
	open($in,'<', $file) or die "$0 read error: file '$file' isn't readable $!";
	}
	$data{'LENMAX'} = 0;
	$data{'QUALMIN'} = 1000;#$data{'QUALMIN'}..$data{'QUALMAX'}
	$data{'QUALMAX'} = 0;
	&ignoreAnythingBeforeFirstFasta(*$in) if($format eq 'FASTA' || $format eq 'MIRDFASTA');
	while(not(eof($in))){
		my $count;
		my $seq;
		my $qual;
		($seq, $count,$qual)= &getNextSequence(*$in, $format);
		chomp($seq);
		($count = 1)if(not($count));
		$seq = reverse($seq)if($reverse);
		$qual = reverse($qual)if($qual && $reverse);
		print STDOUT '#'.$seq.'#'."\n" if($debug);
		$data{'LEN'}{(length($seq)-1)}+=$count;
		($data{'LENMAX'} = length($seq))if($data{'LENMAX'} < length($seq));
		my $counter = 0;
		while($counter < length($seq)){
			my $nuc = uc(substr($seq, $counter,1));
			die "non [ATCGNU] base found in input '$nuc'\n" if(not($nuc=~/[ATCGNU]/));
			$data{'BASE'}{$counter}{$nuc}+=$count;
			$nucleotides{$nuc}++;
			$counter++;
		}
		if($qual){
			if(not(length($seq) == length($qual))){
				print STDERR "Length of qualityscores and sequencelength doesn't correspond:\n$seq\n$qual\n";
				exit(0);
			}
			my $counter = 0;
			while($counter < length($qual)){
				my $qualscore = ord(substr($qual, $counter,1));
				$data{'QUALMIN'} = $qualscore if ($qualscore < $data{'QUALMIN'});
				$data{'QUALMAX'} = $qualscore if ($qualscore > $data{'QUALMAX'});
				$data{'QUAL'}{$counter}{$qualscore}++;
				$counter++;
			}
		}
	}

	close($in);
	print $out $file."\n";
##### printout header
	my $header = "length\tcount\tbase\t".join("\t",sort(keys(%nucleotides)));
	my $header2 = "\t".join("\%\t",sort(keys(%nucleotides)))."\%";
	if($Options{Q}){
		for my $qual($data{'QUALMIN'}..$data{'QUALMAX'}){
			$header2 = $header2."\t".$qual;
		}
	}
	if($Options{q}){
	$header2 = $header2."\t0.10pctqual\t0.50pctqual\t0.90pctqual";
	}
	$header2 = $header2."\n";
	print $out $header."\tsum".$header2;
##### calculate percentages from basecounts and printout both
	for my $base(0..($data{'LENMAX'}-1)){
		print $out ($base+1)."\t";
		print $out $data{'LEN'}{($base)} if($data{'LEN'}{$base});
		print $out "\t".($base+1);
		foreach my $nuc(sort(keys(%nucleotides))){
			print $out "\t";
			if($data{'BASE'}{$base}{$nuc}){
				print $out ($data{'BASE'}{$base}{$nuc});
			}
		}
		
		my $sum;
#		= &arraySum(@data{'BASE'}{$base}{sort(keys(%nucleotides))});
		foreach my $nuc(sort(keys(%nucleotides))){
			($sum += $data{'BASE'}{$base}{$nuc})if($data{'BASE'}{$base}{$nuc});
		}
		print $out "\t".$sum;
		
		foreach my $nuc(sort(keys(%nucleotides))){
			print $out "\t";
			if($sum && $data{'BASE'}{$base}{$nuc}){
				print $out sprintf("%.3f", ($data{'BASE'}{$base}{$nuc})/$sum);#sprintf("%.3f", $number)
			}
		}
##print out quality scores
		#full Qual print;
		if($Options{Q}){
			for my $qual($data{'QUALMIN'}..$data{'QUALMAX'}){
				print $out "\t";
				print $out $data{'QUAL'}{$base}{$qual}if($data{'QUAL'}{$base}{$qual});
			}
		}
		#percentile Qual print
		if($Options{q}){
			my $sumQualCountValsPerBase;
			
			for my $qual($data{'QUALMIN'}..$data{'QUALMAX'}){
				($sumQualCountValsPerBase += $data{'QUAL'}{$base}{$qual})if($data{'QUAL'}{$base}{$qual});
			}
			for my $pct(("0.10","0.50","0.90")){
				#percentile 0.50 mean
				my $pctbasedintervalnum = int($sumQualCountValsPerBase*$pct);
				my $currentintervalnum;
				for my $qual($data{'QUALMIN'}..$data{'QUALMAX'}){
					if($data{'QUAL'}{$base}{$qual}){
						$currentintervalnum += $data{'QUAL'}{$base}{$qual};
					}
					if($currentintervalnum && $pctbasedintervalnum <= $currentintervalnum){
						print $out "\t";
						print $out "$qual";#\t$pctbasedintervalnum\t$currentintervalnum\t$sumQualCountValsPerBase";###check this ... 
					}
					last if($currentintervalnum && $pctbasedintervalnum <= $currentintervalnum);
				}
			}
		}
		#average
		#future work
		
		print $out "\n";
	}
	print STDOUT 'done' if($debug);
}
exit(0);


###############subs and lots off it
# check if input file is fasta, tabdelim, list or fastq format
sub determineFileFormat{
	foreach my $file (@_){
		open(my $handle,'<', $file)or die "$0 read error $file isn't readable $!\n $use";
		my $type = '';
		my $count = 0;
		while(not($type)){
			my $line = <$handle>;
			$type = &checkFasta($line);
			(print STDERR "fas\n")if($debug);
			return $type if($type);
			$type = &checkFastq($line);
			return $type if($type);
			$type = &checkList($line);
			return $type if($type);
			&checkNextline($line, $count);
			return $type if($type);
			$count++;
		}
	}
}
sub ignoreAnythingBeforeFirstFasta{
 #made cause technals issue with reading fast[a]files in this specific way
	local $/ = '>';
	my $fhandle = $_[0];
	my $fastaline = <$fhandle>;
}
#blessed filehandle input and Filetype
sub getNextSequence{
	my $filehandle = $_[0];
	my $format = $_[1];
	if($format eq 'FASTA'){
		&readFasta(*$filehandle);
	}elsif($format eq 'MIRDFASTA'){
		&readmiRDFasta(*$filehandle);
	}elsif($format eq 'FASTQ'){
		&readFastq(*$filehandle);
	}elsif($format eq 'LIST'){
		&readList(*$filehandle);
	}elsif($format =~ /^FTF\d{6}/){
		&readFlexibleTabFormat(*$filehandle,$format);
	}else{
		die "$0 error: invalid file format '$format'\n $!";
	}
}
sub arraySum{
	
	my $sum = 0;
	foreach(@_){
	my $num = $_;
		if($num){
		$sum += $num;
		}
	}
	return $sum;
}
##########################sub determineFileFormat subcheck
	sub checkFasta{
		if($_[0] =~ /^>.*/g){
			#if(#matches miRDeepcollapsed){
				#return 'FASTAmrd';
			#}else{
			return 'FASTA';
			#}
		}
	}
	sub checkFastq{
		if($_[0] =~ /^@.*/){
			#if(#matches fastqcollapsed){
				#return 'FASTQCol';
			#}else{
			return 'FASTQ';
			#}
		}
	}
		sub checkList{
		if($_[0] =~ /^[atcgnATCGN]{1,}/){
			#if($_[0] =~ /^[atcgnATCGN]{1,}\t\d+/){
			#return 'CLIST';
			#}else{
			return 'LIST';
			#}
		}
	}
	sub checkNextline{
		if(($_[0] =~ /#/)&&($_[1] < 200)){
			return '';
		}else{
			die "$0 file format error: Couldn't match the first 200 lines of the file to supported formats $!";
		}
	}
##########################sub getNextSequence subread

	sub readFasta{
		local $/ = '>';
		my $fhandle = $_[0];
		my $fastaline = <$fhandle>;
		chomp($fastaline);
		#print STDOUT $fastaline."\n";
		my @fasta = split("\n",substr($fastaline,1));
		my $seq;
		my $header = shift(@fasta);
		$seq = join('',@fasta);
		$seq =~ s/\r|\n//g;
		return $seq;
	}
	sub readmiRDFasta{
		local $/ = '>';
		my $fhandle = $_[0];
		my $fastaline = <$fhandle>;
		chomp($fastaline);
		#print STDOUT $fastaline."\n";
		my @fasta = split("\n",substr($fastaline,1));
		my $seq;
		my $header = shift(@fasta);
		my $count;
		(undef, $count) = split('_x',$header);
		$seq = join('',@fasta);
		return ($seq,$count);
	}
	sub readFastq{
		my $fileh = $_[0];
		my $header = <$fileh>;
		my $seq = <$fileh>;
		chomp $seq;
		my $header1 = <$fileh>;
		my $qual = <$fileh>;
		chomp $qual;
		return($seq,1,$qual);
	}
	sub readList{
		my $line = <$_[0]>;
		$line =~ /^([atcgnuATCGNU]{1,})/;
		my $seq = $1;
		return $seq;
	}
	sub readFlexibleTabFormat{
		my $line = <$_[0]>;
		my @lineSplit = split('	',$line);
		my $formatKey = $_[1];#format
		if($formatKey =~ /FTF(\d\d)(\d\d)(\d\d)/){
			my $seq = '';
			my $count = 1;
			my $qual = '';
			if(int($1) && int($1) ne 99){
				my $tabseq = int($1);
				$seq = $lineSplit[$tabseq];
			}
			if(int($2) && int($2) ne 99){
				my $tabcount = int($2);
				$count = $lineSplit[$tabcount];
			}
			if(int($3) && int($3) ne 99){
				my $tabqual = int($3);
				$qual = $lineSplit[$tabqual];
			}
			return ($seq,$count,$qual);
		}else{
			die "$0 error:invalid FTF specifier!\n $!";
		}
	}


########################

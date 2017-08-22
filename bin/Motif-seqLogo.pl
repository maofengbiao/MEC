#!/usr/bin/perl -w
#############################################
#Program name:motif-seqLogo
#Author: Fengbiao Mao
#Email:maofengbiao08@163.com || 524573104@qq.com || maofengbiao@gmail.com
#Tel:18810276383
#version:1.1
#Time:Sun Oct 14 22:19:19 2012
##############################################
use strict;
use File::Basename;
use Getopt::Long;
use Cwd qw(abs_path);
use Data::Dumper;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path);
#use PerlIO::gzip;
#use Math::CDF qw(:all);
#use Statistics::Basic qw(:all);
#use Math::Trig;
#############################################
my ($in,$ic,$R,$out,$frag,$motif,$str,$name,$end,$win_max,$Help);
GetOptions
(
	"in=s" => \$in,
	"out=s" => \$out,
	"f=i" => \$frag,
	"c=s" => \$ic,
	"R=s" => \$R,
	"m=s" => \$motif,
	"s=i"=> \$str,
	"e=i"=> \$end,
	"help" => \$Help,
);
##############################################
my $usage=<<USAGE;
\n Usage : \n perl $0  
	-i <enzymed_file>  
	-f <flanking length> [4]
	-m <motif> [CG]
	-c <ic: "T" or "F"> [F]
	-s <start length> [40]
	-e <end length>	[1000]
	-o <out_directory>
	-R <R path> [R]
	-h <display this help info>\n
USAGE
die $usage if ($Help || ! $in || !$out);
#print STDERR "---Program	$0	starts --> ".localtime()."\n";
&showLog("Start!");
##############################################
$str||=40;
$end||=1000;
$frag||=4;
$frag=int($frag);
$motif||="CG";
$ic||="F";
$R||="R";
my $Rlib = "$Bin/../lib/R-packages";
die "Error: -s (start length of fragment) must be not less than 1 && not bigger than '-e' value\n" if ($str<1 || $str>$end);
die "Error: -f (flanking length) must be not less than 1\n" if ($frag<1);
##############################################
if($in=~/\S+.gz$/){
	open IN,"gzip -dc $in |" or die $!;
}else{	
	open IN,"$in" or die $!;
}
if (defined $out)  { open(STDOUT,  '>', "$out/motif_seqLogo.xls")  || die $!; }
my %hash;
#my $all=0;
while(<IN>){
	chomp;
	my @tmp=split;
	next if ($tmp[3] < $str);
	next if ($tmp[3] > $end);
#	next if ($tmp[3] < $frag);
	my $sp=($tmp[3]+1)/$frag;
	while ($tmp[4]=~/(\w{$frag}$motif\w{$frag})/g){
		my $seq=$1;#GCCCCGCCCC
		my @seqs=split ("",$seq);
#		$all++;
		my $n=1;
		foreach my $s (@seqs){
			next if ($s=~/N/i);
			$hash{$n}{$s}++;
			$n++;
		}		
	}		
}
print STDOUT " \tA\tC\tG\tT\n";
foreach my $n (sort {$a<=>$b} keys %hash){
	my @tmp=();
	my $all=0;
	foreach my $s (sort keys %{$hash{$n}}){
		$all+=$hash{$n}{$s};
    }	
	for my $i ("A","C","G","T"){
	#	print STDERR "$i\n";
		if (! exists $hash{$n}{$i}){
			$hash{$n}{$i}=0;
		}
	}
	foreach my $s (sort keys %{$hash{$n}}){
		#print STDERR "$s\n";
		if (! exists $hash{$n}{$s}){
			$hash{$n}{$s}=0;
		}
		my $per=0;
		if($all>0){
			$per=$hash{$n}{$s}/$all;
		}
		push @tmp,$per;
	}
	print STDOUT join "\t",$n,@tmp,"\n";
}
open RS,">$out/Motif.$str-$end.R" or die $!;
print RS "png(\"$out/Motif.$str-$end.png\",height=600,width=800)
rt=read.table(\"$out/motif_seqLogo.xls\",head=T)
library(seqLogo,lib.loc=\"$Rlib\");
seqLogo(t(rt),ic.scale=$ic);
dev.off()
";
`$R CMD BATCH $out/Motif.$str-$end.R`;
&showLog("Done!");

close IN;
sub cal {
    my $winseq=shift @_;
    my $c =$winseq=~s/C/C/g;
    my $g =$winseq=~s/G/G/g;
    my $gc=$winseq=~s/CG/CG/g;
    my $len=length $winseq;
 #   my $gc_frequency=$gc/$len;
    my $gc_ratio=($c+$g)/$len;
        my $gc_OE;
    $gc_OE=$gc*$len/$c/$g;
    return ($gc_ratio,$gc_OE);
}
sub showLog {
        my ($info) = @_;
        my @times = localtime; # sec, min, hour, day, month, year
	print STDERR sprintf("[%d-%02d-%02d %02d:%02d:%02d] %s\n", $times[5] + 1900,$times[4] + 1, $times[3], $times[2], $times[1], $times[0], $info);
}
sub Max{
        my (@aa) = @_;
        my $max=shift @aa;
        foreach  (@aa) {$max=$_ if($_>$max);}
        return $max;
}
#print STDERR "<--Program	$0	ends --- ".localtime()."\n";

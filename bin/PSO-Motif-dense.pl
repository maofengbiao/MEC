#!/usr/bin/perl -w
#############################################
#Program name:PSO-Motif-dense
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
my ($in,$ic,$R,$out,$frag,$motif,$str,$end,$win_max,$Help);
GetOptions
(
	"in=s" => \$in,
	"out=s" => \$out,
	"c=s" => \$ic,
	"R=s" => \$R,
	"m=s" => \$motif,
	"s=f"=> \$str,
	"e=f"=> \$end,
	"help" => \$Help,
);
##############################################
my $usage=<<USAGE;
\n Usage : \n perl $0  
	-i <enzymed_file>  
	-m <motif> [CG]
	-s <start length> [40]
	-e <end length>	[1000]
	-c <horizon or not,T:true,F:False> [F]
	-o <out_directory>
	-R <R path> [R]
	-h <display this help info>\n
USAGE
die $usage if ($Help || ! $in ||!$out);
#print STDERR "---Program	$0	starts --> ".localtime()."\n";
&showLog("Start!");
##############################################
$str||=40;
$end||=1000;
$motif||="CG";
$R||="R";
$ic||="F";
my $Rlib = "$Bin/../lib/R-packages";
die "Error: -s (start length of fragment) must be not less than 1 && not bigger than '-e' value\n" if ($str<1 || $str>$end);
##############################################
if($in=~/\S+.gz$/){
	open IN,"gzip -dc $in |" or die $!;
}else{	
	open IN,"$in" or die $!;
}
if (defined $out)  { open(STDOUT,  '>', "$out/motif_dense_PSO.xls")  || die $!; }
my %hash;
#print STDOUT " \tLength\n";
my $max=0;
while(<IN>){
	chomp;
	my $all=0;
	my @tmp=split;
	next if ($tmp[3] <=0);
	next if ($tmp[3] < $str);
	next if ($tmp[3] > $end);
	while ($tmp[4]=~/($motif)/g){
		$all+=length($motif);
	}
	my $freq=sprintf("%.1f",$all/$tmp[3]);
	if($freq>$max){
		$max=$freq;
	}
	print STDOUT "$freq\t$tmp[3]\n";
}
my @data;
my @name;
for (my $i=0;$i<=$max;$i+=0.1){
	$i=sprintf("%.1f",$i);
	my $d="rt[rt[,1]==$i,2]";
	my $n="\"$i\"";
	push @data,$d;
	push @name,$n;
}
my $dat=join ",",@data;
my $nam=join ",",@name;
#print STDERR"$dat\n$nam\n";
open RS,">$out/Motif-Dense.$str-$end.R" or die $!;
print RS "png(\"$out/Dense.$str-$end.png\",height=600,width=800)
rt=read.table(\"$out/motif_dense_PSO.xls\",head=T)
library(\"sm\",lib.loc=\"$Rlib\");
library(\"vioplot\",lib.loc=\"$Rlib\");
vioplot($dat,rectCol=\"cyan\",names=c($nam),horizontal=$ic)
dev.off()
";
`$R CMD BATCH $out/Motif-Dense.$str-$end.R`;
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

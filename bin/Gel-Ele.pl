#!/usr/bin/perl -w
#############################################
#Program name:Gel-Ele-Fugure
#Author: Fengbiao Mao
#Email:maofengbiao08@163.com || 524573104@qq.com || maofengbiao@gmail.com
#Tel:18810276383
#version:1.1
#Time:Wed Dec 18 23:01:05 2013
##############################################
use strict;
use File::Basename;
use Getopt::Long;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path);
use Data::Dumper;
#use PerlIO::gzip;
#use Math::CDF qw(:all);
#use Statistics::Basic qw(:all);
#use Math::Trig;
#############################################
my ($in,$tran,$str,$end,$R,$name,$min,$max,$out,$cpgpos,$p,$win_max,$Help);
GetOptions
(
	"in=s" => \$in,
	"out=s" => \$out,
	"s=i" => \$str,
	"e=i" => \$end, 
	"t=i" => \$tran,
	"n=s" => \$name,
	"min=i" => \$min, 
	"max=i" => \$max,
	"R=s" => \$R,
	"help" => \$Help,
);
##############################################
my $usage=<<USAGE;
\n Usage : \n perl $0  
	-i <in_file,outfile of Fragment-motifSeq-count.pl>  
	-o <out_directory>
	-s <start length of fragment> [40]
	-e <end length of fragment> [1000]
	-n <names of enzyme> [MspI]
	-t <transparency,0~100> [50]
	-R <R path> [R] 
	-h <display this help info>\n
USAGE
die $usage if ($Help || ! $in || !$out);
&showLog("Start!");
##############################################
my @files=`ls $in`;
my $fh;
$str||=40;
$tran||=50;
$end||=1000;
$R||="R";
$name||="MspI";
die "Error: -s (start length of fragment) must be not less than 1 && not bigger than '-e' value\n" if ($str<1 || $str>$end);
die "Error: -t (transparency) must be not less than 0\n" if ($tran <0);
##############################################
if (defined $out)  { open(STDOUT,  '>', "$out/Gel.$str-$end.R")  || die $!; }
foreach my $file (@files){
	chomp $file;
	print STDOUT "png(\"$out/Gel.$str-$end.png\")
par(mar = c(1, 1, 1, 1))
par(bg = \"grey30\")
bin=$tran
plot(20, log(3000), type = \"n\", ylim = c(log(30), log(3000)),bty='n',
xlim = c(0, 10), yaxt = \"n\", xaxt = \"n\")
a = c(10,50,seq(100, 800, by = 100), 1500, 3000)
colMar = paste(\"grey\", c(10,20,35, 50, 55, 60, 70, 100, 80, 85,90, 100), sep = \"\")
text(1.5, log(a), paste(a, \"bp\",sep=\"\"), col = \"white\")
segments(2.5 - 0.28, log(a), 2.5 + 0.28, log(a), col = colMar,lwd = 4)
rt=read.table(\"$file\")
#first
a=rt[rt[,1]>=$str & rt[,1]<=$end,1]
colst = ceiling(100*rt[rt[,1]>=$str & rt[,1]<=$end,2]/max(rt[rt[,1]>$str & rt[,1]<=$end,2]))
cols =ceiling(100*(colst+bin)/(100+bin))
col1 = paste(\"grey\", cols, sep = \"\")
i=2
segments((1.5 + i) - 0.28, log(a), (1.5 + i) + 0.28, log(a), col = col1, lwd = 10)
text(3.5, log(3000), \"Motif\", col = \"white\",font = 2,cex=1.2)
#second
colst = ceiling(100*rt[rt[,1]>$str&rt[,1]<=$end,3]/max(rt[rt[,1]>$str&rt[,1]<=$end,3]))
cols =ceiling(100*(colst+bin)/(100+bin))
col1 = paste(\"grey\", cols, sep = \"\")
i=3
segments((1.5 + i) - 0.28, log(a), (1.5 + i) + 0.28, log(a), col = col1, lwd = 10)
text(4.5, log(3000), \"DNA\", col = \"white\",font = 2,cex=1.2)
#third
rt1=rt[rt[,1]>$str&rt[,2]>0&rt[,1]<=$end,2]/rt[rt[,1]>$str&rt[,2]>0&rt[,1]<=$end,3]
colst = ceiling(100*rt1/max(rt1))
cols =ceiling(100*(colst+bin)/(100+bin))
col1 = paste(\"grey\", cols, sep = \"\")
i=4
segments((1.5 + i) - 0.28, log(a), (1.5 + i) + 0.28, log(a), col = col1, lwd = 10)
text(5.5, log(3000), \"M/D\", col = \"white\",font = 2,cex=1.2)
text(7, log(500), \"$name\", col = \"white\",font = 2,cex=1.5)
";
}
`$R CMD BATCH $out/Gel.$str-$end.R`;
&showLog("Done!");
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
sub min{
        my (@aa) = @_;
        my $min=shift @aa;
        foreach  (@aa) {$min=$_ if($_<$min);}
        return $min;
}
sub open_file {
        my $file=shift;
		my $fh;
        if($file=~/\S+.gz$/){
                open $fh,"gzip -dc $file |" or die $!;
        }else{
        open $fh,"$file" or die $!;
        }
        return $fh;
}
#print STDERR "<--Program	$0	ends --- ".localtime()."\n";

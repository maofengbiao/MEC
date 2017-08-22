#!/usr/bin/perl -w
#############################################
#Program name:fragment_motifSeq_count
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
#use PerlIO::gzip;
#use Math::CDF qw(:all);
#use Statistics::Basic qw(:all);
#use Math::Trig;
#############################################
my ($in,$out,$R,$motif,$str,$name,$end,$win_max,$Help);
GetOptions
(
	"in=s" => \$in,
	"out=s" => \$out,
	"n=s" => \$name,
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
	-m <motif> [CG]
	-n <name of enzymes> [MspI]
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
$motif||="CG";
$str||=40;
$end||=1000;
$name||="MspI";
die "Error: -s (start length of fragment) must be not less than 1 && not bigger than '-e' value\n" if ($str<1 || $str>$end);
##############################################
if($in=~/\S+.gz$/){
	open IN,"gzip -dc $in |" or die $!;
}else{	
	open IN,"$in" or die $!;
}
if (defined $out)  { open(STDOUT,  '>', "$out/fragment_count.xls")  || die $!; }
my %hash;
while(<IN>){
	chomp;
	my @tmp=split;
	my $len=length($motif);
	my $cg=$len*($tmp[4]=~s/$motif/$motif/g);
	$hash{$tmp[3]}[0]+=$cg;
	$hash{$tmp[3]}[1]+=$tmp[3];
}
foreach my $key (sort {$a<=>$b} keys %hash){
	print  STDOUT "$key\t$hash{$key}[0]\t$hash{$key}[1]\n";
}
open RS,">$out/Frag.$str-$end.R" or die $!;
print RS "png(\"$out/Frag.$str-$end.png\",height=600,width=800)
a<-read.table(\"$out/fragment_count.xls\")
x=a[,1]
y=a[,3]/(sum(as.numeric(a[,3])))
par(font.lab=2,font.axis=1,cex.lab=1,cex.axis=1,mar=c(3.5,4.2,2,3.5),mgp=c(2,0.7,0))
plot(x,100*y,ylim=c(0,max(100*y)),xlim=c($str,$end),main=\"$name\",ylab=\"\",xlab=\"Fragment Length\",type=\"l\",pch=19,col=\"green3\",lwd=2)
mtext(\"Percentage of Fragment Divide by Genome(%)\",side = 2,line = 2,cex=1,font=2,col=\"green3\")
par(new=T)
x=a[a[,3]>0,1]
y=a[a[,3]>0,2]/a[a[,3]>0,3]
plot(x,100*y,ylim=c(0,max(100*y)),xlim=c($str,$end),ylab=\"\",xlab=\"\",axes=F,type=\"l\",pch=19,col=\"BlueViolet\",lwd=2)
axis(4,cex=2,cex.lab=2,font=1) 
mtext(\"Mean Density of $motif in Fragments of Same Length (%)\",side = 4,line = 2,cex=1,font=2,col=\"BlueViolet\")
";
`$R CMD BATCH $out/Frag.$str-$end.R`;
&showLog("Done!");

close IN;
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

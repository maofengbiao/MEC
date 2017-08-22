#!/usr/bin/perl -w
#############################################
#Program name:PSO
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
my ($in,$str,$end,$R,$min,$max,$out,$cpgpos,$p,$win_max,$Help);
GetOptions
(
	"in=s" => \$in,
	"out=s" => \$out,
	"s=i" => \$str,
	"e=i" => \$end, 
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
	-min <min distance between start and end> [200]
	-max <min distance between start and end> [500]
	-R <R path> [R] 
	-h <display this help info>\n
USAGE
die $usage if ($Help || ! $in || !$out);
#print STDERR "---Program	$0	starts --> ".localtime()."\n";
&showLog("Start!");
##############################################
my @files=`ls $in`;
my $fh;
$str||=40;
$end||=1000;
$min||=200;
$max||=500;
$R||="R";
die "Error: -s (start length of fragment) must be not less than 1 && not bigger than '-e' value\n" if ($str<1 || $str>$end);
die "Error: -min (min distance between start and end) must be not less than 1 && not bigger than '-max' value\n" if ($min<1 || $min>$max);
##############################################
if (defined $out)  { open(STDOUT,  '>', "$out/pso.$str-$end.R")  || die $!; }
my $Rlib = "$Bin/../lib/R-packages";
foreach my $file (@files){
	chomp $file;
	print STDOUT "beginLength=$str
endLength=$end
minWide=$min
maxWide=$max
library(pso,lib.loc=\"$Rlib\")
set.seed(1)
setwd(\"$out\")
data<-read.table(\"$file\",sep=\"\\t\",header=FALSE)
data<-sapply(data,as.numeric)
targetFunc<-function(x){
        x1=as.integer(x[1])
        x2=x1+x[2]

        if(x1>=x2)
                return(.Machine\$integer.max) 
        if((x2-x1<minWide)||(x2-x1>maxWide))
                        return(.Machine\$integer.max)    
        lineBegin=which(data[,1]>=x1)[1]
        lineEnd=which(data[,1]>=x2)[1]
        return(-1*sum(data[lineBegin:lineEnd,2])/sum(as.numeric(data[lineBegin:lineEnd,3])))
}
myPso<-psoptim(rep(NA,2),targetFunc,lower=c(beginLength,minWide),upper=c(endLength,maxWide),control=list(maxit=1000,s=10))
flag1t=myPso\$par[1]
flag2t=myPso\$par[2]+myPso\$par[1]
flag1=paste(flag1t,\"bp\",sep=\"\")
flag2=paste(flag2t,\"bp\",sep=\"\")
write.table(data.frame(flag1,flag2),file=\"$out/pso.result.xls\",sep=\"\\t\",quote = F,col.names=c(\"Start\",\"End\"),row.names=F)
#write.table(paste(flag1,flag2,sep=\"\\t\"),file=\"$out/pso.result.xls\",sep=\"\\t\");
"
}
`$R CMD BATCH $out/pso.$str-$end.R`;
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

#!/usr/bin/perl -w
#############################################
#Program name:GC_OE_count
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
die "Error: -s (start length of fragment) must be not less than 1 && not bigger than '-e' value\n" if ($str<1 || $str>$end);
##############################################
if($in=~/\S+.gz$/){
	open IN,"gzip -dc $in |" or die $!;
}else{	
	open IN,"$in" or die $!;
}
if (defined $out)  { open(STDOUT,  '>', "$out/gc_oe_count.xls")  || die $!; }
my %hash;
while(<IN>){
	chomp;
	my @tmp=split;
	next if ($tmp[3] < $str);
	next if ($tmp[3] > $end);
	my ($gc,$oe)=cal($tmp[4]);	
	$gc=sprintf("%.1f",$gc);
	$oe=sprintf("%.1f",$oe);
	$hash{$gc}[0]++;
	$hash{$oe}[1]++;
}
foreach my $key (sort {$a<=>$b} keys %hash){
	if (! exists $hash{$key}[0]){$hash{$key}[0]=0;}
	if (! exists $hash{$key}[1]){$hash{$key}[1]=0;}
	print  STDOUT "$key\t$hash{$key}[0]\t$hash{$key}[1]\n";
}

open RS,">$out/GC-OE.$str-$end.R" or die $!;
print RS "png(\"$out/GC-OE.$str-$end.png\",height=600,width=800)
rt=read.table(\"$out/gc_oe_count.xls\")
dat<-data.frame(rt[,2]/sum(rt[,2]),rt[,3]/sum(rt[,3]))
row.names(dat)<-sprintf(\"%.1f\",rt[,1])
#row.names(dat)<-seq(0.1,1,by=0.1)
a<-barplot(t(as.matrix(dat)),xlab=\"Content\",las=2,font.lab=2,
ylab=\"Frequency\",
main=\"GC Content and CpGoe Value Distribution\",
legend=c(\"GC-content\",\"CpGoe-value\"),
yaxt=\"n\",col=c(\"green3\",\"BlueViolet\"),
ylim=c(0,1.1*max(dat)),beside=T,
args.legend=list(bty=\"n\",cex=1.2));
sp<-spline(a[1,],dat[,1],n=200)
lines(sp\$x,sp\$y,lwd=2,col=\"green3\")
sp<-spline(a[2,],dat[,2],n=200)
lines(sp\$x,sp\$y,lwd=2,col=\"BlueViolet\");
axis(2,seq(0,1,by=0.1),labels=c(\"0%\",\"10%\",\"20%\",\"30%\",\"40%\",\"50%\",\"60%\",\"70%\",\"80%\",\"90%\",\"100%\"))
";
`$R CMD BATCH $out/GC-OE.$str-$end.R`;
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

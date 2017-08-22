#!/usr/bin/perl -w
#############################################
#Program name:Base-distribution
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
my ($in,$R,$out,$frag,$motif,$str,$name,$end,$win_max,$Help);
GetOptions
(
	"in=s" => \$in,
	"out=s" => \$out,
	"f=i" => \$frag,
	"n=s" => \$name,
	"m=s" => \$motif,
	"s=i"=> \$str,
	"e=i"=> \$end,
	"R=s" => \$R,
	"help" => \$Help,
);
##############################################
my $usage=<<USAGE;
\n Usage : \n perl $0  
	-i <enzymed_file>  
	-f <Relative fragments> [100]
	-s <start length> [40]
	-e <end length>	[1000]
	-R <R path> [R]
	-o <out_directory>
	-h <display this help info>\n
USAGE
die $usage if ($Help || ! $in || !$out);
#print STDERR "---Program	$0	starts --> ".localtime()."\n";
&showLog("Start!");
##############################################
$str||=40;
$end||=1000;
$frag||=100;
$frag=int($frag);
die "Error: -s (start length of fragment) must be not less than 1 && not bigger than '-e' value\n" if ($str<1 || $str>$end);
die "Error: -f (Relative fragments) must be bigger than 1\n" if ($frag<1); 
##############################################
if($in=~/\S+.gz$/){
	open IN,"gzip -dc $in |" or die $!;
}else{	
	open IN,"$in" or die $!;
}
if (defined $out)  { open(STDOUT,  '>', "$out/base_percent_count.xls")  || die $!; }
my %hash;
while(<IN>){
	chomp;
	my @tmp=split;
	next if ($tmp[3] < $str);
	next if ($tmp[3] > $end);
	next if ($tmp[3] < $frag);
	my $sp=($tmp[3]+1)/$frag;
	my $n=1;
	for (my $i=1;$i<$tmp[3]+1;$i+=$sp){
		my $s=int($i);
		my $t=int($i+$sp);
		my $l=$t-int($i);
#		die "length ERROR\n" if ($l<=0);
		my $seq=substr($tmp[4],$s-1,$l);
	#	print STDERR "$tmp[3]\t$sp\t$i\t$s\t$t\t$l\n";
#		die "Sequence ERROR\n" if ($seq eq "");
		my $As=($seq=~s/A/A/ig);
		my $Ts=($seq=~s/T/T/ig);
		my $Cs=($seq=~s/C/C/ig);
		my $Gs=($seq=~s/G/G/ig);
		my $all=$As+$Ts+$Cs+$Gs;
		$hash{$n}[0]+=$As;
		$hash{$n}[1]+=$Ts;
		$hash{$n}[2]+=$Cs;
		$hash{$n}[3]+=$Gs;
		$hash{$n}[4]+=$all;
		$n++;
	}
}
print  STDOUT " \tA\tT\tC\tG\n";
foreach my $key (sort {$a<=>$b} keys %hash){
	if (! exists $hash{$key}[0]){$hash{$key}[0]=0;}
	if (! exists $hash{$key}[1]){$hash{$key}[1]=0;}
    if (! exists $hash{$key}[2]){$hash{$key}[2]=0;}
    if (! exists $hash{$key}[3]){$hash{$key}[3]=0;}
	if (! exists $hash{$key}[4]){$hash{$key}[4]=0;}
	my $Ap=$hash{$key}[0]/$hash{$key}[4];
	my $Tp=$hash{$key}[1]/$hash{$key}[4];
	my $Cp=$hash{$key}[2]/$hash{$key}[4];
	my $Gp=$hash{$key}[3]/$hash{$key}[4];
	print  STDOUT "$key\t$Ap\t$Tp\t$Cp\t$Gp\n";
}

open RS,">$out/Bases.$str-$end.R" or die $!;
print RS "png(\"$out/Bases.distribution.$str-$end.png\",height=600,width=800)
rt=read.table(\"$out/base_percent_count.xls\",head=T)
barplot(t(as.matrix(rt)),ylim=c(0,1.1),xlab=\"Relative Position (5'->3')\",ylab=\"Relative Frequency\",main=\"Bases Distribution\",pch=15,col=c(\"mediumorchid\",\"cyan\",\"lightgreen\",\"sandybrown\"),font.lab=2);
legend(\"top\",legend=c(\"A\",\"T\",\"C\",\"G\"),pch=15,col=c(\"mediumorchid\",\"cyan\",\"lightgreen\",\"sandybrown\"),horiz=T,bty=\"n\",cex=1.4)
";
`$R CMD BATCH $out/Bases.$str-$end.R`;
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

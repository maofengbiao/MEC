#!/usr/bin/perl -w
#############################################
#Program name:Reports
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
my ($in,$l,$pso,$mec,$out,$str,$enz,$name,$end,$win_max,$Help);
GetOptions
(
	"in=s" => \$in,
	"out=s" => \$out,
	"t=s" => \$pso,
	"l=s" => \$l,
	"n=s" => \$enz,
	"m=s" => \$mec,
	"s=i"=> \$str,
	"e=i"=> \$end,
	"help" => \$Help,
);
##############################################
my $usage=<<USAGE;
\n Usage : \n perl $0  
	-i <Figure directory by Magic Enzyme Cutter>
	-t <The table of PSO results>  
	-l <Logo of Icons>
	-m <Logo of MEC>
	-n <names of enzymes> [MspI]
	-s <start length> [40]
	-e <end length>	[1000]
	-o <out_directory>
	-h <display this help info>\n
USAGE
die $usage if ($Help || ! $in || !$out);
#print STDERR "---Program	$0	starts --> ".localtime()."\n";
&showLog("Start!");
##############################################
$str||=40;
$end||=1000;
$enz||="MspI";
die "Error: -s (start length of fragment) must be not less than 1 && not bigger than '-e' value\n" if ($str<1 || $str>$end);
##############################################
if (defined $out)  { open(STDOUT,  '>', "$out/Report-MEC.html")  || die $!; }
my $fh=open_file($pso);
my ($t1,$t2,$t3,$t4);
while(<$fh>){
	chomp;
	my @tmp=split;
	if($.==1){
		($t1,$t2)=(@tmp)[0,1];
	}
	if($.==2){
        ($t3,$t4)=(@tmp)[0,1];
    }
}
close $fh;
my @time = localtime;
my $date=sprintf("[%d-%02d-%02d %02d:%02d:%02d]", $time[5] + 1900,$time[4] + 1, $time[3], $time[2], $time[1], $time[0]);
print STDOUT "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Strict//EN\">
<html>
<head><title>MEC MspI Report</title>

<style type=\"text/css\">

 \@media screen {
  div.summary {
    width: 18em;
    position:fixed;
    top: 10em;
    margin:1em 0 0 1em;
  }
  
  div.main {
    display:block;
    position:absolute;
    overflow:auto;
    height:auto;
    width:auto;
    top:11.8em;
    bottom:2.3em;
    left:18em;
    right:0;
    border-left: 1px solid #CCC;
    padding:0 0 0 1em;
    background-color: white;
    z-index:1;
  }
  
  div.header {
    background-color: #EEE;
    border:0;
    margin:0;
    padding: 0.5em;
    font-size: 200%;
    font-weight: bold;
    position:fixed;
    width:100%;
    top:0;
    left:0;
    z-index:2;
  }

  div.footer {
    background-color: #EEE;
    border:0;
    margin:0;
	padding:0.5em;
    height: 1.3em;
	overflow:hidden;
    font-size: 100%;
    font-weight: bold;
    position:fixed;
    bottom:0;
    width:100%;
    z-index:2;
  }
  
  img.indented {
    margin-left: 3em;
  }
 }
 
 \@media print {
	img {
		max-width:100% !important;
		page-break-inside: avoid;
	}
	h2, h3 {
		page-break-after: avoid;
	}
	div.header {
      background-color: #FFF;
    }
	
 }
 
 body {    
  font-family: sans-serif;   
  color: #000;   
  background-color: #FFF;
  border: 0;
  margin: 0;
  padding: 0;
  }
  
  div.header {
  border:0;
  margin:0;
  padding: 0em;
  font-size: 200%;
  font-weight: bold;
  width:100%;
  }    
  
  #header_title {
  display:inline-block;
  float:left;
  clear:left;
  font-size: 80%;
  margin-top:0em;
  }
  #header_filename {
  display:inline-block;
  float:right;
  clear:right;
  font-size: 50%;
  margin-right:3em;
  text-align: right;
  }

  div.header h3 {
  font-size: 50%;
  margin-bottom: 0;
  }
  
  div.summary ul {
  padding-left:0;
  list-style-type:none;
  }
  
  div.summary ul li img {
  margin-bottom:-0.5em;
  margin-top:0.5em;
  }
	  
  div.main {
  background-color: white;
  }
      
  div.module {
  padding-bottom:1.5em;
  padding-top:1.5em;
  }
	  
  div.footer {
  background-color: #EEE;
  border:0;
  margin:0;
  padding: 0.5em;
  font-size: 100%;
  font-weight: bold;
  width:100%;
  }


  a {
  color: #000080;
  }

  a:hover {
  color: #800000;
  }
      
  h2 {
  color: #800000;
  padding-bottom: 0;
  margin-bottom: 0;
  clear:left;
  }

  table { 
  margin-left: 3em;
  text-align: center;
  }
  
  th { 
  text-align: center;
  background-color: #000080;
  color: #FFF;
  padding: 0.4em;
  }      
  
  td { 
  font-family: monospace; 
  text-align: left;
  background-color: #EEE;
  color: #000;
  padding: 0.4em;
  }

  img {
  padding-top: 0;
  margin-top: 0;
  border-top: 0;
  }

  
  p {
  padding-top: 0;
  margin-top: 0;
  }
  
</style>

</head>

<body>

<div class=\"header\">
<div id=\"header_title\"><img  align=center src=\"$mec\" alt=\"Magic-Enzyme-Cutter\" width=\"170\" height=\"170\"> <font color=darkgreen> Magic-Enzyme-Cutter Report&nbsp </div>
<p><font size=\"3\" font color=dimgray> Master Restriction Enzyme II digestion in silico && Evaluate the degradation products in many aspects && Calculate the optimal selection in expect (MEC).<br />Enzymes: $enz<br />&nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp $date</font></p>
</div>
<div class=\"summary\">
<h2>Summary</h2>
<ul>
<li><img src=\"$l\" alt=\"[PASS]\"> <a href=\"#M0\">Fragment && motifSeq count</a></li>
<li><img src=\"$l\" alt=\"[PASS]\"> <a href=\"#M1\">Gel electrophoresis imaging</a></li>
<li><img src=\"$l\" alt=\"[PASS]\"> <a href=\"#M2\">GC && OE count</a></li>
<li><img src=\"$l\" alt=\"[PASS]\"> <a href=\"#M3\">Bases distribution</a></li>
<li><img src=\"$l\" alt=\"[PASS]\"> <a href=\"#M4\">SeqLogo of Motif</a></li>
<li><img src=\"$l\" alt=\"[PASS]\"> <a href=\"#M5\">PSO  optimization</a></li>
<li><img src=\"$l\" alt=\"[PASS]\"> <a href=\"#M6\">Motif density after PSO</a></li>
<a href=\"Digest/Digested_Genome.txt\"> <br /> Digested-Genome.txt</a> presents a detailed information of each fragment. Format is \"chromosome start end length sequence\".</p>
</div>

<div class=\"main\">
<div class=\"module\"><h2 id=\"M0\"><img src=\"$l\" alt=\"[OK]\"> Fragment && motifSeq count</h2>
<p><img class=\"indented\" src=\"$in/Frag.$str-$end.png\" alt=\"Fragment && motifSeq count graph\" width=\"720\" height=\"560\"></p>
</div>
<div class=\"module\"><h2 id=\"M1\"><img src=\"$l\" alt=\"[OK]\"> Gel electrophoresis imaging</h2>
<p><img class=\"indented\" src=\"$in/Gel.$str-$end.png\" alt=\"Gel electrophoresis imaging graph\" width=\"720\" height=\"560\"></p>
</div>
<div class=\"module\"><h2 id=\"M2\"><img src=\"$l\" alt=\"[OK]\"> GC && OE count</h2>
<p><img class=\"indented\" src=\"$in/GC-OE.$str-$end.png\" alt=\"GC && OE count graph\" width=\"720\" height=\"560\"></p>
</div>
<div class=\"module\"><h2 id=\"M3\"><img src=\"$l\" alt=\"[OK]\"> Bases distribution</h2>
<p><img class=\"indented\" src=\"$in/Bases.distribution.$str-$end.png\" alt=\"Bases distribution graph\" width=\"720\" height=\"560\"></p>
</div>
<div class=\"module\"><h2 id=\"M4\"><img src=\"$l\" alt=\"[OK]\"> SeqLogo of Motif</h2>
<p><img class=\"indented\" src=\"$in/Motif.$str-$end.png\" alt=\"SeqLogo of Motif graph\" width=\"720\" height=\"560\"></p>
</div>
<div class=\"module\"><h2 id=\"M5\"><img src=\"$l\" alt=\"[OK]\"> PSO  optimization</h2>
<table>
<tr>
<th>Start</th>
<th>End</th>
</tr>
<tr>
<td>40bp</td>
<td>240bp</td>
</tr>
</table>
</div>
<div class=\"module\"><h2 id=\"M6\"><img src=\"$l\" alt=\"[OK]\"> Motif density after PSO</h2>
<p><img class=\"indented\" src=\"$in/Dense.40-240.png\" alt=\"Motif density after PSO graph\" width=\"720\" height=\"560\"></p>
</div>
</div>
";
&showLog("Done!");

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

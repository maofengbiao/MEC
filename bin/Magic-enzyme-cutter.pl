#!/usr/bin/perl 
use strict;
use warnings;
use FindBin qw($Bin $Script);
if (@ARGV<3){
    die "\nperl $0 <ref.fa> <type1>,<type2>... <brief code> [ outfile ]
example:
    perl $0 hg19.fa  C-CGG,C-GGC,G-CGC     brief-code.txt  > outfile.txt
    perl $0 hg19.fa  MspI_AciI+TaqI        brief-code.txt  > outfile.txt
    perl $0 hg19.fa  MspI+TTT-AAA+YAC-GTR  brief-code.txt  > outfile.txt\n\n";
}
my $file = shift;
my $type=shift;
my $code=shift;
my $outfile=shift;
if (defined $outfile){
	open STDOUT,">$outfile" or die $!;
}
my $chr;
my $seq;

##############--------begin-------################
open NEB,"$Bin/enzyme.list.good.txt" or die $!;
my %neb;
while(<NEB>){
	chomp;#AatII   GACGT-C
	my @tmp=split;
	$neb{$tmp[0]}=$tmp[1];
}
my @types=split (/[\,\_\+]/,$type);
my %enzymes;
my $flag=0;
foreach my $t (@types){
	chomp $t;
	my $tp;
	my $pos;
	if($t!~/\-/){#name of NEB
		if(exists $neb{$t}){
			$tp=$neb{$t};	
			$pos = index($neb{$t},"-"); 
		}else{
			die "\nError: no such restriction endonuclease name $t, please check the spelling (Case Sensitivity) or use motif sequence such as \"C-CGG\"\n\n";
		}
	}else{#motif sequence
		$tp=$t;
		$pos = index($t,"-"); 
	}
	$flag++;
	$tp=~s/\-//;
	$tp=uc($tp);
	$enzymes{$flag}=[$pos,$tp];#pos=0 ; # typeseq=1
}
if($file=~/\S+.gz/){
        open IN,"gzip -dc $file |" or die $!;
}else{  
        open IN,"$file" or die $!;
}

if($code=~/\S+.gz/){
        open CD,"gzip -dc $code |" or die $!;
}else{
        open CD,"$code" or die $!;
}
#--code of NEB---#
my %hash;#neb coding
while(<CD>){#R=G/A
	chomp;
	next if ($_=~/^\s*$/);
	next if ($_=~/^#/);
	my $line=$_;
	if($line=~/(\w)=(\S+)/){
		my $name=$1;
		my @tmp=split(/\//,$2);
		$hash{$name}=\@tmp;
	}
}
my %patterns;
foreach my $f (sort {$a<=>$b} keys %enzymes){
	my @abcd=split //,$enzymes{$f}[1];
	my @all;
	foreach my $abc (@abcd){
		if ($abc!~/[ATCG]/){#brief code
			if ($hash{$abc}){
				push @all,[@{$hash{$abc}}];
			}else{
				die "\nError: do not exist such code for RE from NEB ...$abc\n\n";
			}
		}else{
			push @all,[$abc];
		}
	}
	my $list=\@all;
	my $str_list=create($list);
	my @pat=@$str_list;	#for one enyzme;
	$patterns{$f}=\@pat;
}
#foreach my $it (@{$str_list}){
#        print STDERR "$it\tenzyme-seq\n";
#}
#store reference
$/ = ">"; <IN>; $/ = "\n";
while (<IN>){
	chomp;
	$chr = (split ( /\s+/, $_))[0];
#	print STDERR "chr-->$chr\n";
	$/ = ">";
	$seq = <IN>;
	$/ = "\n";
	chomp($seq);
	$seq =~ s/\s+//g;
	$seq =~ s/>$//;
	$seq = uc($seq);
	deal ($chr,\$seq);
}
close (IN);

##################----sub program----#################
sub create {
    my ($list) = @_;
    my $str_list = [ '' ];
    foreach my $ref (@{$list}) {
        $str_list = create_list($str_list, $ref);
    }
    return $str_list;
}

sub create_list {
    my ($ref_str_list, $ref_array) = @_;
    my @return_array;
    foreach my $str (@{$ref_str_list}) {
        foreach my $element (@{$ref_array}) {
            push @return_array, "${str}$element";
        }
    }
    return \@return_array;
}
##########
sub deal{
	my ($chr,$seq) = @_;
	my %allstart1;
	my %allstart2;
	my $start1;
	my $start2;
	#get the min start1 of every enzymes
	foreach my $f (sort {$a<=>$b} keys %enzymes){
	    my @start1=(); 
		my @patt=@{$patterns{$f}};
		foreach my $ele (@patt){#foreach each patten
	#		print STDERR "patten-->$ele\n";
			my $st1 = index($$seq,"$ele");
	#		print STDERR "start1-->$st1\n";
			push @start1,$st1;
		}
		my $minstart1 = min(@start1);
		$allstart1{$f}=$minstart1;
	}
	#get the min start1 from all enzymes
    my $tag=0;
    my $flag1=0;
    foreach my $k (keys %allstart1){
        if ($tag==0){#frist one
            $start1=$allstart1{$k};
            $flag1 = $k; 
            $tag=1;
        }else{#next one
            if ($allstart1{$k}<$start1){
                $start1 = $allstart1{$k};
                $flag1 = $k;    
            }   
        }   
    }   
	#get the min start2 of every enzymes
	foreach my $f (sort {$a<=>$b} keys %enzymes){	
		my @start2=();
		my @patt=@{$patterns{$f}};
		foreach my $ele (@patt){ 
            my $st2 = index ($$seq,"$ele",$start1 + 1);
            push @start2,$st2;
        }
		my $minstart2 = min(@start2);
		$allstart2{$f}=$minstart2;
	}
	if ($start1 == -1){#no any match; exit to out
		my $chrlength=length($$seq);
		print STDOUT "$chr\t1\t$chrlength\t$chrlength\t$$seq\n";
		return;
	}
	#get the min start2 from all enzymes
	my $tag2=0;
	my $flag2=0;
	foreach my $k (keys %allstart2){
        if ($tag2==0){
            $start2=$allstart2{$k};
			$flag2 = $k;
            $tag2=1;
        }else{
			if ($allstart2{$k}<$start2){
				$start2 = $allstart2{$k};
				$flag2 = $k;
			}
        	}
    	} 
	my $length;
	my $seqs;
	my $pos=$enzymes{$flag1}[0];
	my $tmp_len=$start1 + $pos;
	$seqs = substr ($$seq,0,$tmp_len);
	print STDOUT join "\t",$chr,1,$tmp_len,$tmp_len,$seqs;
	print STDOUT "\n";
	while ($start2 != -1){ #correction
#		my @start2=();
		$length = $start2 - $start1;
		$seqs = substr ($$seq,$start1 + $pos,$length);
		print STDOUT join "\t",$chr,$start1+$pos+1,$start2 + $pos,$length,$seqs;
		print STDOUT "\n";
		$start1 = $start2;
		### new run ###
		foreach my $f (sort {$a<=>$b} keys %enzymes){
        	my @start2=();
       		my @patt=@{$patterns{$f}};
        	foreach my $ele (@patt){
            	my $st2 = index ($$seq,"$ele",$start1 + 1);
            	push @start2,$st2;
        	}
        	my $minstart2 = min(@start2);
        	$allstart2{$f}=$minstart2;
    	}	
	my $tag2=0;
	my $flag2=0;
    	foreach my $k (keys %allstart2){
        	if ($tag2==0){
            	$start2=$allstart2{$k};
            	$flag2 = $k;
            	$tag2=1;
        	}else{
            	if ($allstart2{$k}<$start2){
                	$start2 = $allstart2{$k};
                	$flag2 = $k;
            	}
        	}
    	}
		$pos=$enzymes{$flag2}[0]; 
    }
	$seqs = substr ($$seq,$start1 + $pos);
	my $len_final=length($seqs);
	my $tmp_end=$len_final+($start1 + $pos);
	print STDOUT join "\t",$chr,$start1 + $pos+1,$tmp_end,$len_final,$seqs;
	print STDOUT "\n";
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

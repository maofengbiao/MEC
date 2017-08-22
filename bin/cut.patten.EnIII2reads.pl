#!/usr/bin/perl 
use strict;
use warnings;
my $file = shift;
my $type=shift;
my $code=shift;
my $expand=shift;
my $read=shift;
my $site=shift;
my $outfile=shift;
my $chr;
my $seq;
die "perl $0 <ref.fa> <type>(C-CGG) <code> <expand length> <readout> <sites file> [ outfile ]\n" unless ($file);
open READ,">$read" or die $!;
open SITE,">$site" or die $!;
if (defined $outfile){
    open STDOUT,">$outfile" or die $!;
}
my $tp=$type;
my $pos = index($type,"-");
print STDERR "pos of index cut site\t$pos\n";
$tp=~s/\-//;
$tp=uc($tp);
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
my %hash;
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
my @abcd=split //,$tp;
my @all;
foreach my $abc (@abcd){
	if ($abc!~/[ATCG]/){
		if ($hash{$abc}){
			push @all,[@{$hash{$abc}}];
		}else{
			die "Error, do not exist such code for RE from NEB ...$abc\n";
		}
	}else{
		push @all,[$abc];
	}
}
#foreach my $k (@all){
#	print STDERR "@$k\n";
#}
my $list=\@all;
my $str_list=create($list);
my @patterns=@$str_list;

foreach my $it (@{$str_list}){
        print STDERR "$it\tenzyme-seq\n";
}
$/ = ">"; <IN>; $/ = "\n";
while (<IN>){
	chomp;
	$chr = (split ( /\s+/, $_))[0];
	$/ = ">";
	$seq = <IN>;
	$/ = "\n";
	chomp($seq);
	$seq =~ s/\s+//g;
	$seq =~ s/>$//;
	$seq = uc($seq);
	&deal ($chr,\$seq,\@patterns,$expand);
}
close (IN);

#########################################
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
sub deal {
	my ($chr,$seq,$pat,$expand) = @_;
	my @patterns=@$pat;
	my @start1=();
	my @start2=();
	print STDERR length($$seq)."\n";
	foreach my $t (@patterns){
		my $st1 = index ($$seq,"$t");
		push @start1,$st1;
	}
	my $start1 = min(@start1);#the first motif site
	my $start11 = $start1+$expand;
	my $start1b= $start1+1-$expand;
	if($start1b<1){$start1b=1;}
	foreach my $t (@patterns){
        my $st2 = index ($$seq,"$t",$start11 + $pos);
        push @start2,$st2;
    }
	my $start2 = min(@start2);#the second motif site
	my $start22 = $start2+$expand;
	my $start2b = $start2+1-$expand;
	if($start2b<1){$start2b=1;}
	my $length;
	my $seqs;
	my $start11bk;
	my $tmp_len=$start1b-1;
	if($tmp_len>0){
		$seqs = substr ($$seq,0,$tmp_len);
		print STDOUT join "\t",$chr,1,$tmp_len,$tmp_len,$seqs,"\n";#begin of cutting
	}
	my $len=($start11+$pos)-($start1b)+1;
	$seqs = substr ($$seq,$start1b-1,$len);

	my $len1=($start1+1-1+2)-($start1b)+1;
    my $read1= substr ($$seq,$start1b-1,$len1);
	if ($len1>0){
    	print READ ">Left-$chr-".$start1b.":".($start1+1-1+2)."@".$chr."-".$start1b.":".($start11+$pos)."\n";
        print READ "$read1\n";
    }
    my $len2=($start11+$pos)-($start1+$pos+1-2)+1;
    my $read2=substr ($$seq,$start1+$pos+1-2-1,$len2);
    if ($len2>0){
        print READ ">Right-$chr-".($start1+$pos+1-2).":".($start11+$pos)."@".$chr."-".$start1b.":".($start11+$pos)."\n";
    	print READ "$read2\n";
    }
	print STDOUT join "\t",$chr,$start1b,$start11+$pos,$len,$seqs,"\n";
	print SITE join "\t",$chr,$start1+1+1,"\n";
	while ($start2 != -1){#motif2
		my @start2=();
		if ($start2b<=$start11+$pos){
			$start2b=$start11+$pos+1;
		}else{ #internal of two motifs
			$len=$start2b-1-($start11+$pos+1)+1;
			if($len>0){
				$seqs = substr ($$seq,$start11+$pos+1-1,$len);
				print STDOUT join "\t",$chr,$start11+$pos+1,$start2b-1,$len,$seqs,"\n";
			}
		}
		$length=($start22+$pos)-($start2b)+1;
		$seqs = substr ($$seq,$start2b-1,$length);

		my $len1=($start2+1-1+2)-($start2b)+1;
		my $read1= substr ($$seq,$start2b-1,$len1);
		if ($len1>0){
			print READ ">Left-$chr-".$start2b.":".($start2+1-1+2)."@".$chr."-".$start2b.":".($start22+$pos)."\n";
			print READ "$read1\n";
		}
		my $len2=($start22+$pos)-($start2+$pos+1-2)+1;
		my $read2=substr ($$seq,$start2+$pos+1-2-1,$len2);
		if ($len2>0){
			print READ ">Right-$chr-".($start2+$pos+1-2).":".($start22+$pos)."@".$chr."-".$start2b.":".($start22+$pos)."\n";
			print READ "$read2\n";
		}
		print STDOUT join "\t",$chr,$start2b,$start22+$pos,$length,$seqs,"\n";
		print SITE join "\t",$chr,$start2+1+1,"\n";
#iteration
		$start11bk=$start11;
		$start11 = $start22;
		foreach my $t (@patterns){
        	my $st2 = index ($$seq,"$t",$start11 + $pos );
        	push @start2,$st2;
    	}
		$start2 = min(@start2);
		if($start2 != -1){
			$start22 = $start2+$expand;	
			$start2b = $start2+1-$expand;
		}
    }
# last fragment
	$seqs = substr ($$seq,$start22+$pos);#substr to the end
	my $len_final=length($seqs);
	if ($len_final>0){
		my $tmp_end=$len_final+($start22+$pos);
		print STDOUT join "\t",$chr,$start22+$pos+1,$tmp_end,$len_final,$seqs,"\n";
	}
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

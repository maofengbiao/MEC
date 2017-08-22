-------------------------
Magic-Enzyme-Cutter (MEC)
-------------------------

Author : 
	Mao FengBiao,Cai WanShi,Wen YanLing.
Contact: 
	maofengbiao@gmail.com.
Version: 
	1.0.
License:
	GPLv3 (http://www.gnu.org/licenses/gpl.html)

--------------------
Software description
--------------------

	Magic-Enzyme-Cutter is a software that Master Restriction Enzyme II
digestion in silico && Evaluate the degradation products in many aspects &&  
Calculate the optimal selection in expect (MEC).

------------------------
Platform and environment
------------------------
    
 System: Linux or UNIX
 Software: perl >= v5.8.8, R >= R-2.15.1 are required.

perl download website : http://www.perl.org/
R download website : http://www.r-project.org/
 
--------------
How to run MEC
--------------

1. Install the software

>tar xzvf MEC-xxx.tar.gz
>cd MEC-xxx
>sh install /path/to/R

2. Run the software

>cd MEC-xxx
>./Magic-Enzyme-Cutter -help

3. Test the software

>cd MEC-xxx/test
>sh test.sh /path/to/R

Results of test are in directory 'MEC-xxx/test/test-result'

More information is available in file 'MEC-xxx/Manual' and 
in directory 'MEC-xxx/doc' which contains *doc files or PDF 
files about introduction and manual.

-------------------
Magic-Enzyme-Cutter
-------------------

        Master Restriction Enzyme II digestion in silico && Evaluate the
degradation products in many aspects && Calculate the optimal selection in
expect (MEC).
        Author : Mao FengBiao,Cai WanShi,Wen YanLing.
        Contact: maofengbiao@gmail.com.
        Version: 1.0.
Usage:

        Magic-Enzyme-Cutter -i <ref.fa> -o <outdir> [options] 

Require:
        -i <in_file: reference in FASTA format>
        -o <out_directory>

Options:
        -R <R path> [R]
        -z <enzyme name or cut-Type> [MspI]
        -m <motif which is focus on> [CG]
        -s <start length of fragment> [40]
        -e <end length of fragment> [1000]
        -n <min regions of PSO search> [200]
        -x <max regions of PSO search> [500]
        -f <bases of flanking included> [4]
        -g <graduated scale number> [100]
        -c <change the vision,T or F> [F]
        -t <transparency of Gel Figure> [50]
        -h <display this help information>

Note:
        The delimiter of enzymes is ",", "_" or "+" while the delimiter of
cut-Type is "-".
        The brief-codes of cut-Type are supported, They are:
R=G/A;Y=C/T;M=A/C;K=G/T;S=G/C;W=A/T;B=C/G/T;D=A/G/T;H=A/C/T;V=A/C/G;N=A/C/G/T.
        Gzipped FASTA format input is supported.

Example:
        Magic-Enzyme-Cutter -i <ref.fa> -o <outdir> -z MspI -s 100 -e 2000
        Magic-Enzyme-Cutter -i <ref.fa> -o <outdir> -R /Path/R -z C-CGG,C-GGC,G-CGC
        Magic-Enzyme-Cutter -i <ref.fa> -o <outdir> -R /Path/R -z MspI_AciI+TaqI
        Magic-Enzyme-Cutter -i <ref.fa> -o <outdir> -R /Path/R -z MspI+TTT-AAA+YAC-GTR
        Magic-Enzyme-Cutter -i <ref.fa> -o <outdir> -n 100 -x 300 -g 50 -f 5 -c T 

The FASTA format file of sequences. See example data:

>chrR
GATCTGATAAGTCCCAGGACTTCAGAAGagctgtgagaccttggccaagt
cacttcctccttcagGAACATTGCAGTGGGCCTAAGTGCCTCCTCTCGGG
ACTGGTATGGGGACGGTCATGCAATCTGGACAACATTCACCTTTAAAAGT
TTATTGATCTTTTGTGACATGCACGTGGGTTCCCAGTAGCAAGAAACTAA
AGGGTCGCAGGCCGGTTTCTGCTAATTTCTTTAATTCCAAGACAGTCTCA
AATATTTTCTTATTAACTTCCTGGAGGGAGGCTTATCATTCTCTCTTTTG
GATGATTCTAAGTACCAGCTAAAATACAGCTATCATTCATTTTCCTTGAT
TTGGGAGCCTAATTTCTTTAATTTAGTATGCAAGAAAACCAATTTGGAAA
TATCAACTGTTTTGGAAACCTTAGACCTAGGTCATCCTTAGTAAGATctt
cccatttatataaatacttgcaagtagtagtgccataattaccaaacata
aagccaactgagatgcccaaagggggccactctccttgcttttcctcctt

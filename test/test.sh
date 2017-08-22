#!/bin/bash
if [ $# -lt 1 ]
then
        echo ""
        echo " Usage: $0 <R path>"
        echo ""
        exit
fi
R=$1;
pwd_dir=`pwd`
perl $pwd_dir/../Magic-Enzyme-Cutter  -i  ref.fa.gz  -o  $pwd_dir/test-result -R $R

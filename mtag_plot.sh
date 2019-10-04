#!/bin/bash


###
###
###	requires /usr/local/bin/regional_plot.R
###	plink must be called by 'plink2' 
###	

time=$(date)
echo ""
echo "### filename: mtag_plot - generate manhattan, qq, and regional plots from mtag output"
echo "### author: Matthew Traylor"
echo "### email: m.traylor@qmul.ac.uk"
echo "### last updated 13/10/2019"
echo ""
echo "Analysis started at $time"

if [[ $# -lt 5 ]] ; then
    echo ''
    echo 'Insufficient arguments entered. Please enter the following six arguments:'
    echo '1/ input file '
    echo '2/ output filename prefix'
    echo '3/ p-value threshold for local plots to be drawn'
    echo '4/ location of reference data (prefix only)'
    echo '5/ location of gene-range file (full file location)'
    echo '6/ location of genetic map files (location only)'
    echo 'exiting...'
    echo ''
    exit 1
fi

awk 'OFS="\t"{print $1,$4,$5,$8,$6,$9,$10,$12,$11,$7,$2,$3,"NA","NA","NA"}' $1 > $2 

echo "generating QQ, manhattan and QC plots..."
R --vanilla --args $2 < ./manhattan_qq.R > analysis.out

echo "extracting regions with p<'$3'..."
awk '{if(NR==1)print "SNP","P";else if($8<'f')print $1,$8}' f=$3 $2 > $2.sig_regions


echo "clumping on significant regions..."
plink2 --silent --bfile $4 --clump $2.sig_regions --clump-r2 0.1 --clump-kb 1000 --out $2.sig_regions >> analysis.out 


echo "extracting significant regions and plotting..."
sed '/^$/d' $2.sig_regions.clumped > tmp && mv tmp $2.sig_regions.clumped 
ln=$(wc -l $2.sig_regions.clumped | awk '{print $1}')
i=2
while [[ $i -le $ln ]]
do
		j=$((i-1))
        SNP=$(awk 'NR=='f'{print $3}' f=$i $2.sig_regions.clumped)
        BP=$(awk 'NR=='f'{print $4}' f=$i $2.sig_regions.clumped)
        CHR=$(awk 'NR=='f'{print $1}' f=$i $2.sig_regions.clumped)

        let "BP_min=$BP-250000"
        let "BP_max=$BP+250000"

	plink2 --bfile $4 --ld-snp $SNP --r2 --ld-window-kb 250 --ld-window-r2 0 --ld-window 99999 --out $2.sig_region"$j" >> analysis.out
        awk '{if($11='f' && $12>'g' && $12<'h')print $1,$11,$12,$2,$3,$4,$5,$6,$7,$8}' f=$CHR g=$BP_min h=$BP_max $2 > $2.sig_region"$j"
        ((i=i+1))
        
        R --vanilla --args $2.sig_region"$j" $CHR $2.sig_region"$j".ld $5 $6< ./regional_plot.R >> analysis.out

done
time=$(date)
echo "Analysis completed at $time"

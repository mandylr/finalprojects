#! usr/bin/env bash

#This script will take a DMR and methylation bed file and make a generate
#a bigwig for use in UCSC browswer

#Select only DMR's with a fdr < .01

#-log($4)/log(10) is the same as -log10($4)

cat data/DMR.bed | cut -f2,3,4,5,6,7 \
    | awk '$4 < .00000001' \
    | awk 'BEGIN{OFS="\t"} $4 = -log($4)/log(10)'\
    | cut -f1,2,3,4 \
    | bedtools sort -i - > results/DMR_sig.bed

bedGraphToBigWig results/DMR_sig.bed ../data-sets/genome/hg19.genome \
    results/DMR.bw

#Make bedgraph of meth.bed for CBD and BeS values

cat data/meth.bed | tail -n +2 |cut -f1,2,3,5 \
    | bedtools sort -i - > results/meth_BeS.bg
    
bedGraphToBigWig results/meth_BeS.bg ../data-sets/genome/hg19.genome \
    results/meth_BeS.bw

cat data/meth.bed | tail -n +2 | cut -f1,2,3,6 \
    | bedtools sort -i - > results/meth_CBD.bg

bedGraphToBigWig results/meth_CBD.bg ../data-sets/genome/hg19.genome \
    results/meth_CBD.bw

#Get genes assocaiated with DMRs
#Downloaded genes.bed in RefGene reference from UCSC genome table Genes
#and Gene Prediction

bedtools intersect -a results/DMR_sig.bed -b ../genes.sort.bed -wb \
    | awk -v OFS="\t" '{print $1,$2,$3,$4,$8,$9}' \
    > results/DMR.genes.bed

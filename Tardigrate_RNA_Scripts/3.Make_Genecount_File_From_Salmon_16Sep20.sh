############

cut -f1 T_1_quant.sf > genename
for f in T_1_quant.sf T_2_quant.sf T_3_quant.sf T_4_quant.sf T_5_quant.sf T_6_quant.sf T_7_quant.sf T_8_quant.sf T_9_quant.sf T_10_quant.sf T_11_quant.sf T_12_quant.sf T_13_quant.sf T_14_quant.sf T_15_quant.sf T_16_quant.sf T_17_quant.sf T_18_quant.sf 
do
    sample=$f
    awk 'BEGIN {OFS="\t"} {print $5}' $f | sed -e "1s/NumReads/$sample/" > $f.tmp
done
paste genename T_1_quant.sf.tmp T_2_quant.sf.tmp T_3_quant.sf.tmp T_4_quant.sf.tmp T_5_quant.sf.tmp T_6_quant.sf.tmp T_7_quant.sf.tmp T_8_quant.sf.tmp T_9_quant.sf.tmp T_10_quant.sf.tmp T_11_quant.sf.tmp T_12_quant.sf.tmp T_13_quant.sf.tmp T_14_quant.sf.tmp T_15_quant.sf.tmp T_16_quant.sf.tmp T_17_quant.sf.tmp T_18_quant.sf.tmp > genecount

###########
# convert a file to .csv
cat genecount | tr -s '[:blank:]' ',' > genecount.csv

############

cut -f1 T_1_quant.sf > genename
for f in T_1_quant.sf T_2_quant.sf
do
    sample=$f
    awk 'BEGIN {OFS="\t"} {print $5}' $f | sed -e "1s/NumReads/$sample/" > $f.tmp
done
paste genename T_1_quant.sf.tmp T_2_quant.sf.tmp > genecount

###########
# convert a file to .csv
cat genecount | tr -s '[:blank:]' ',' > genecount.csv

#! /bin/bash

echo -e "Desc\tFamily ID\t$(head -1 HOG_stat_filtered.tsv |cut -f2-)" > gene_families_cafe.txt

col=$(head -1 HOG_stat_filtered.tsv|tr '\t' '\n'|grep -n -w "DROSOPHILA_MELANOGASTER"|cut -d ":" -f1)

cat ../no_lambda_hog.txt > tmp.txt
#cat control_hog.txt >> tmp.txt
# Orthologs with >=10 species and must contain genes from Dmel
# Root level hog only
awk -v FS='\t' -v OFS='\t' '{non_empty=0; for(i=2; i<=NF; i++) if($i >= "1") non_empty++; if(non_empty >= 5) print $0}' <(tail -n +2 HOG_stat_filtered.tsv)|\
awk -v FS='\t' -v OFS='\t' -v a=125 '($a > 0) {print $0}'| awk -v FS='\t' -v OFS='\t' '{
    min = $2
    max = $2
    for(i=3; i<=NF; i++) {
        if($i < min) min = $i
        if($i > max) max = $i
    }
    if(max - min < 100) print $0
}'| awk '{print "(null)\t" $0}' >> gene_families_cafe.txt

#! /bin/bash

# This script retrieves the number of interactions (union of physical and genetic) for each HOG in the given file.

# Usage: ./get_interactions_perHOG.sh control_genes.tsv

> interactions_perHOG.tsv
 echo -e "HOG\tgroup\tgenetic\tphysical" > interactions_perHOG.tsv

>HOG_not_found.tsv

get_interactions_perHOG() {
    og=$1
    FBgn=$2
    group=$3   
genetic=$(sort -u genetic_interactions_DroID_flybase.tsv|awk -v gene="$FBgn" '{
    if ($1 == gene && $2 != gene) arr[$2]=1;
    if ($2 == gene && $1 != gene) arr[$1]=1;
}
END {
    n=0;
    for (g in arr) {
        printf n++ ? ", %s" : "%s", g
    }
    print ""
}' )


physical=$(sort -u protein_protein_interactions_DroID.tsv|awk -v gene="$FBgn" '{
    if ($1 == gene && $2 != gene) arr[$2]=1;
    if ($2 == gene && $1 != gene) arr[$1]=1;
}
END {
    n=0;
    for (g in arr) {
        printf n++ ? ", %s" : "%s", g
    }
    print ""
}' )

    echo -e "${og}\t${group}\t${genetic}\t${physical}" >> interactions_perHOG.tsv
    
}

export -f get_interactions_perHOG

while read -r I_gene hog C_genes C_trans; do

FBgn=$(grep -w "gene-Dmel_"${I_gene} NP_FBid.tsv|cut -f3|sort -u)

get_interactions_perHOG "$hog" "$FBgn" "$hog"

# Convert comma-separated strings to arrays
IFS=',' read -ra arr1 <<< "$C_genes"
IFS=',' read -ra arr2 <<< "$C_trans"

# Loop over the arrays in parallel
for i in "${!arr1[@]}"; do
        C_tran=${arr2[$i]}
        C_gene=$(grep -w "${C_tran}" NP_FBid.tsv | cut -f2|sed 's/gene-Dmel_//')
        C_hog=$(grep -w ${C_tran} ../HOG_filtered.tsv|cut -f1)
        C_FBgn=$(grep -w "${C_tran}" NP_FBid.tsv | cut -f3)

        # hog is empty 
        if [ -z "$C_hog" ]; then
            echo -e "${C_gen}\t${C_tran}" >> HOG_not_found.tsv
        else
        get_interactions_perHOG "$C_hog" "$C_FBgn" "$hog"
        fi
done

done < <(tail -n +2 "$1"| cut -f1-3,5)
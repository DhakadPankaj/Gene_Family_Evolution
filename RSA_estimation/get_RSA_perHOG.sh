#! /bin/bash


# get RSA for each HOG, and proportions of exposed and buried residues (RSA > 0.25 and RSA < 0.1, respectively)
# Usage: ./get_RSA_perHOG.sh control_genes.tsv
> RSA_perHOG.tsv
echo -e "HOG\tRSA\tExposed\tBuried\tgroup" > RSA_perHOG.tsv

get_rsa() {
    hog=$2
    I_gene=$1
    C_genes=$3
    C_trans=$4

get_RSA_perHOG() {
    og=$1
    gn=$2
    group=$3

    # Get average RSA value for the gene
    read RSA exposed buried < <(cat rsa_values/${gn}.rsa.csv|tr "," "\t"|cut -f5|tail -n +2|awk '{sum+=$1; n+=1; if($1>=0.25) p+=1; if($1<=0.1) s+=1} END {print sum/n, p/n, s/n}')
    
    if [ -z "$RSA" ]; then
        RSA="NA"
        exposed=0
        buried=0
    fi

    echo -e "${og}\t${RSA}\t${exposed}\t${buried}\t${group}" >> RSA_perHOG.tsv
}

export -f get_RSA_perHOG


#while read -r I_gene hog C_genes C_trans; do

    FBgn=$(grep -w "gene-Dmel_"${I_gene} NP_FBid.tsv | cut -f3|sort -u)

    pdb=$(grep -w "${I_gene}" uniprotkb_proteome_available.tsv)
    if [ -z "$pdb" ]; then
        echo "No PDB entry found for $I_gene"
    else
    if [ "$(echo "$pdb" | wc -l)" -gt 1 ]; then
        pdb=$(echo "$pdb" | awk -F'\t' 'NR==1{first=$0} $2=="reviewed"{print; found=1; exit} END{if(!found) print first}')
    fi

    pdb=$(echo "$pdb" | cut -f1)

    if [ -s "pdb/AF-${pdb}-F1-model_v4.pdb" ]; then
        python3 calc_rsa.py pdb/AF-${pdb}-F1-model_v4.pdb -o rsa_values/${I_gene}
        get_RSA_perHOG "$hog" "$I_gene" "$hog"
        else
        echo -e "${hog}\tNA\tNA\tNA\tNA" >> RSA_perHOG.tsv
    fi

    fi
    # Convert comma-separated strings to arrays
    IFS=',' read -ra arr1 <<< "$C_genes"
    IFS=',' read -ra arr2 <<< "$C_trans"

    # Loop over the arrays in parallel
    for i in "${!arr1[@]}"; do
        C_tran=${arr2[$i]}
        C_gene=$(grep -w "${C_tran}" NP_FBid.tsv | cut -f2|sed 's/gene-Dmel_//')
        C_hog=$(grep -w ${C_tran} ../HOG_filtered.tsv|cut -f1)

        C_FBgn=$(grep -w "${C_tran}" NP_FBid.tsv | cut -f3)

        C_pdb=$(grep -w "${C_gene}" uniprotkb_proteome_available.tsv)
        if [ -z "$C_pdb" ]; then
            echo "No PDB entry found for $C_gene"
        else
            if [ "$(echo "$C_pdb" | wc -l)" -gt 1 ]; then
                C_pdb=$(echo "$C_pdb" | awk -F'\t' 'NR==1{first=$0} $2=="reviewed"{print; found=1; exit} END{if(!found) print first}')
            fi

            C_pdb=$(echo "$C_pdb" | cut -f1)

            if [ -s "pdb/AF-${C_pdb}-F1-model_v4.pdb" ]; then
            
                python3 calc_rsa.py pdb/AF-${C_pdb}-F1-model_v4.pdb -o rsa_values/${C_gene}
                get_RSA_perHOG "$C_hog" "$C_gene" "$hog"
                else
                echo -e "${C_hog}\tNA\tNA\tNA" >> RSA_perHOG.tsv
            fi
        fi
    done

#done < <(tail -n +2 "$1"| cut -f1-3,5)
}

export -f get_rsa
parallel -j 1 --colsep "\t" get_rsa :::: <(tail -n +2 "$1"| cut -f1-3,5)
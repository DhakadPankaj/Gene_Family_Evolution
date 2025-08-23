#! /bin/bash

# Run HyPhy analysis on gene families (MEME and FUBAR)

run_hyphy() {

OGid=$1
if [ -s og_alignments/${OGid}/${OGid}_final_align_NT.aln ] && [ -s Gene_trees/${OGid}.treefile ]; then

for method in MEME FEL; do
mkdir -p HyPhy_${method}
echo "HYPHYMPI ${method} --alignment og_alignments/${OGid}/${OGid}_final_align_NT.aln --tree Gene_trees/${OGid}.treefile --output HyPhy_${method}/${OGid}_${method}.json"
HYPHYMPI ${method} --alignment og_alignments/${OGid}/${OGid}_final_align_NT.aln --tree Gene_trees/${OGid}.treefile CPU=10 --output HyPhy_${method}/${OGid}_${method}.json > HyPhy_${method}/${OGid}_output.txt 2> HyPhy_${method}/${OGid}_error.txt
done
fi
}

export -f run_hyphy
parallel -j 6 --env run_hyphy --colsep "\t" run_hyphy {1} :::: hog_hyphy_MEME.txt

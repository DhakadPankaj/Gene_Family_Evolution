#! /bin/bash

#$ -N HyPhy
#$ -cwd
#$ -V
#$ -t 1003-2284
#$ -tc 100
#$ -pe smp64 4
#$ -e std_out/
#$ -o std_out/

source /ceph/users/pdhakad/miniconda3/etc/profile.d/conda.sh && conda activate hyphy_env

method="BUSTED"

OGid=`cat cafe_gene_family.txt|tail -n +2|cut -f2 | awk "NR==$SGE_TASK_ID"`

mkdir -p HyPhy_${method}

# Redirect stdout and stderr to the desired files
exec > HyPhy_${method}/${OGid}_output.txt 2> HyPhy_${method}_error.txt

if [ -s og_alignments/${OGid}/${OGid}_final_align_NT.aln ] && [ -s GeneTrees/${OGid}.treefile ]; then
echo "HYPHYMPI ${method} --alignment og_alignments/${OGid}/${OGid}_final_align_NT.aln --tree GeneTrees/${OGid}.treefile --output HyPhy_${method}/${OGid}_${method}.json"
HYPHYMPI ${method} --alignment og_alignments/${OGid}/${OGid}_final_align_NT.aln --tree GeneTrees/${OGid}.treefile CPU=4 --output HyPhy_${method}/${OGid}_${method}.json
fi

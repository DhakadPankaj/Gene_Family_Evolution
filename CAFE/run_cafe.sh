#! /bin/bash

basename -a family_files/*tsv | parallel -j 6 'cafe5 -i family_files/{} -t hog_UCLDtree.nw -c 10 -p -k 2  --lambda_per_family -eBase_error_model.txt -o failed_cafe/{.} > failed_cafe/{.}.out'

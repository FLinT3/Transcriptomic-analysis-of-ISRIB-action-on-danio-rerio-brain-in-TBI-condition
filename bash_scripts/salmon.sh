c=1

while [ $c -le 27 ]; do

mkdir QUANTS_1/${c}

salmon quant \
-i indexes \
-l A \
-p 16 \
--validateMappings \
--gcBias \
 --numGibbsSamples 23 \
 -o QUANTS_1/${c} \
 -1 SRR_R1_001.fastq.gz \
 -2 SRR_R2_001.fastq.gz

((c ++ ))

done

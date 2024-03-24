trimmomatic PE \
-threads 24 \
input_forward.fq.gz \ 
input_reverse.fq.gz \
output_forward_paired.fq.gz \
output_forward_unpaired.fq.gz \
output_reverse_paired.fq.gz \
output_reverse_unpaired.fq.gz \
SLIDINGWINDOW:2:20 \
MINLEN:60 \
TRAILING:5 \
ILLUMINACLIP:Illumina_PE.fa:2:32:16

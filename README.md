# AIPDS
conda create -n blast

conda activate blast

conda install bioconda::blast

blastp -query new_sequences.fasta -db sabdab_cdr_seq_label_all_db -outfmt 6 -out blast_results.txt

or 

python calculate_blast_similarity.py #change the input

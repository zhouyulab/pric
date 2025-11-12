alt_genome=$1
PREFIX=$2
Gh_genome="../Supplymental/Ghir_genome/Ghirsutum_HAU_genome.fasta"
# step1 
python scripts/extract_site_flank_fasta.py -g $Gh_genome -f 500 -o results/$2_case_flank.fasta

# step2
minimap2 -a ../Supplymental/Ghir_genome/TM_1.min results/$2_case_flank.fasta > results/$2_case_flank.sam
samtools view -F 3844 -bS -o results/$2_case_flank.bam results/$2_case_flank.sam

### step3
python scripts/alignment_pos.py -f results/$2_case_flank.fasta \
-g $Gh_genome \
-b results/$2_case_flank.bam \
-o results/$2_case_pair_position.txt -s results/$2_case_pair_sequence.txt

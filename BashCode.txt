CODE Github repository (BeetGenomeNinja)
Trimmomatic
java -jar $TRIM/trimmomatic PE ${File}_R1_001.fastq.gz ${File}_R2_001.fastq.gz ${File}_R1_pe.fastq ${File}_R1_se.fastq ${File}_R2_pe.fastq ${File}_R2_se.fastq ILLUMINACLIP:/mnt/scratch/galewski/CROP/CT_EL_SR_PAT_genomes/admera_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36




Alignment and VCF:
ssh galewski@hpc.msu.edu
/mnt/home/GATK/Readme
bwa mem -t 10 /mnt/j/EL10_genome_files/EL10byChromosome.fasta XXX XXX > XXX.bam
gatk AddOrReplaceReadGroups -I KDH13.bam -O KDH13.out.bam -LB String -PL String -PU String -SM String
RGadd.sh
gatk AddOrReplaceReadGroups -I KDH13.bam -O KDH13.out.bam -LB String -PL String -PU String -SM String
SamSort.sh
gatk SortSam –I .out.bam -O out.sorted.bam -SO coordinate
BamIndex.sh
gatk BuildBamIndex -I KDH13.out.bam -O KDH13.out.bam.bai
GATK_HapC.sh
gatk HaplotypeCaller -R ../EL10_genome_files/EL10byChromosome.fasta -I KDH13.out.bam -O KDH13.out.vcf
Filter Genotypes (GVCF)
Bash filters ?
Analysis Ready Variants (GVCF) indels and SNP variation
MERGE VCFFILES
F1,P1
CTs
KDH13
Plot unique variation onto genome






BLASTP:
tblastn -db Beta_vulgaris.fa -query KDH13.protiens.fa -num_threads 4 -outfmt 7 > KDH13_to_EL10_2.out
cat KDH13_PacBio.pep | grep gene | grep AUG > peptde_KDH13_order
Where do these proteins blast into the genome?
cat KDH13_PacBio.pep | grep gene | grep '30_259314415_contig_1       ' | awk '{print $9}' | sed 's/ID=//g'
 
cat KDH13toEL10.blast | grep -a1 hit | grep -v -e '--' | grep EL | grep -v '#' > KDH13toEL10_1_besthits.blast
 
for i in $(awk '{print $1}' contig_genes.txt | uniq);do grep $i'     ' contig_genes.txt | tail -1 >> contig_genes1.txt;done
 
augustus KDH13_PacBio.fa --species=arabidopsis --strand=both --gff3=on > KDH13_PacBio.pep
 
cat KDH13_PacBio.pep | grep '#' | grep -v end | grep -v '###' | sed 's/# //g' | sed 's/protein sequence = \[//g' | sed 's/]//g'
 
cat KDH13_PacBio.pep | grep '#' | grep -v end | grep -v '###' | sed 's/# //g' | sed 's/protein sequence = \[//g' | sed 's/]//g' | grep -v '#' | grep -v '-'
 
cat KDH13_PacBio.pep | grep '#' | grep -v end | grep -v '###' | sed 's/# //g' | sed 's/protein sequence = \[//g' | sed 's/]//g' | grep -v '#' | grep -v '-' | grep -v Pred
 
vcftools --relatedness2 --vcf KDH_F1_P1.vcf --out KDH_F1_P1.realtedness
 
(py3) paul@arsidkim4dPJG9K:/mnt/z/2020/BLAST$ cat KDH13toEL10.blast | awk '$3 > "75"' | grep -v '#' | awk '{print $1}' | uniq | wc -l
27847 genes
 
(base) paul@arsidkim4dPJG9K:/mnt/z/2020/BLAST/KDH13pep_slpit$ awk '{print $1}' pfam_scan_out_sorted.txt | uniq | wc -l
40779 abinitio genes with predicted protein domains (Pfam)
 
(py3) paul@arsidkim4dPJG9K:/mnt/z/2020/BLAST$ cat KDH13toEL10.blast | awk '$3 > "75"' | grep -v '#' | awk '{print $2}' | sort | uniq | wc -l
18009
 
View VCF KPF1
awk '{print $1,$2,$3,$4,$5,$6,$9,$10,$11,$12}' KDH_F1_P1.vcf | grep 0/1 | grep 1/1 | cgrep -l A,T,G,C
 
 
cat KDH_F1_P1.vcf | awk '{print $10}' | grep 0/1 | wc -l
2098218 variants were het, genotyping error, when depth was enough and GQ was high, potential
Double haploid suggesting problems in alignment across 20% of the genome. Repeats in beets.
 
 
module load bcftools
bcftools convert -o <output.vcf.gz> -O z <input.vcf>
module load tabix
tabix -p vcf L19.all.vcf.gz
bcftools merge CROPTYPE_FILES > CROPTYPE.vcf
for i in $(ls);do bcftools convert -o $i.gz -O z $i;done

#gunzip -c data/Scerevisiae.fasta.gz >Scerevisiae.fasta
#samtools faidx Scerevisiae.fasta
#samtools dict Scerevisiae.fasta >Scerevisiae.dict
#bwa index Scerevisiae.fasta

#wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR101/SRR10178655/SRR10178655.sra

#fastq-dump --gzip --split-files --defline-seq '@$sn/$ri' --defline-qual '+' SRR10178655.sra

#fastqc --expand --outdir . SRR10178655_1.fastq.gz SRR10178655_2.fastq.gz

#java -Xmx500m -jar ./software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -summary SRR10178655.summary SRR10178655_1.fastq.gz SRR10178655_2.fastq.gz SRR10178655_1_passed.fastq SRR10178655_1_failed.fastq SRR10178655_2_passed.fastq SRR10178655_2_failed.fastq MINLEN:100 ILLUMINACLIP:./data/adapters.fa:2:30:10:2:keepBothReads

bwa mem -t 4 -R '@RG\tID:SRR10178655\tSM:SRR10178655\tLB:SRR10178655\tPL:ILLUMINA' ./Scerevisiae.fasta SRR10178655_1_passed.fastq SRR10178655_2_passed.fastq | samtools view -u - | samtools sort -m 1G -o SRR10178655.bam -

samtools index SRR10178655.bam

gatk MarkDuplicates --java-options '-Xmx1G'  --VALIDATION_STRINGENCY SILENT -I SRR10178655.bam -O SRR10178655.mdup.bam -M SRR10178655.mdup.metrics -MAX_FILE_HANDLES 2000

samtools stats --ref-seq Scerevisiae.fasta SRR10178655.mdup.bam >SRR10178655.stats
plot-bamstats -s Scerevisiae.fasta >Scerevisiae.gc
plot-bamstats -r Scerevisiae.gc -p SRR10178655 SRR10178655.stats

samtools view -b -f3 -F3852 SRR10178655.mdup.bam > SRR10178655.mdup.filter.bam
samtools index SRR10178655.mdup.filter.bam

gatk HaplotypeCaller --minimum-mapping-quality 30 --min-base-quality-score 20 --read-validation-stringency SILENT --reference Scerevisiae.fasta --input SRR10178655.mdup.filter.bam --output SRR10178655.vcf

# using the 'hist' script that we wrote...
samtools depth -q 20 -Q 30 SRR10178655.mdup.filter.bam | hist 3 1 - >SRR10178655.mdup.filter.depth.hist

vcftools --recode --recode-INFO-all --stdout --max-missing 1 --minDP 21 --maxDP 57 --vcf SRR10178655.vcf >SRR10178655.filter.vcf

zcat ~/pfb2019/workshops/NGS/data/Scerevisiae.gff3.gz | grep -E '^#|CDS' >cds.gff3

bedtools intersect -nonamecheck -u -b cds.gff3 -a SRR10178655.filter.vcf | wc -l  # 6180

grep -v '^#' SRR10178655.filter.vcf | wc -l  # 10625


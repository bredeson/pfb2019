NGS Workshop Tutorial
=====================

**NOTE:** Unless otherwise specified, example command lines are available in this workshop's [lecture notes](https://github.com/bredeson/pfb2019/blob/master/workshops/NGS/bio_info_formats.pdf) 

1. First, install the `gnuplot` command line plotting software using HomeBrew:
   ```bash
   brew install gnuplot
   ```

2. Download the sofware required for this tutorial using `wget` from `https://www.dropbox.com/s/zn4m0oqbasg4gal/software.tar.gz`. Unarchive the files by typing `tar -xzvf software.tar.gz`. Finally, add the `software/bin` directory to your `$PATH` environment variable:
   ```bash
   export PATH=$PWD/software/bin:$PATH;
   ```
   **NOTE**: You need to replaced `$PWD` above with the full, expanded path to the `software/bin` dir if you put the `export` statement in your `.bash_profile` or `.bashrc`.
   
   **ANSWER**: Download the software package with `wget`:
   ```bash
   wget https://www.dropbox.com/s/zn4m0oqbasg4gal/software.tar.gz
   tar -xzvf software.tar.gz
   ```
   You can add the `software` directory to your `.bash_profile` or `.bashrc` like so (assuming I downloaded `software.tar.gz` to `/Users/jessenbredeson/pfb2019/workshops/NGS`:
   ```bash
   export PATH=/Users/jessenbredeson/pfb2019/workshops/NGS:$PATH
   ```

3. Download the [genome](https://www.dropbox.com/s/goo2bt4br9mqxtt/Scerevisiae.fasta.gz) and [annotation](https://www.dropbox.com/s/uq8mfp125jlgknq/Scerevisiae.gff3.gz) files using `wget` and decompress them both with `gunzip`.

4. Index the genome as described in the lecture notes.
   
   **ANSWER**: Be sure to check that your downloaded file is *actually* a FASTA file by looking at it with `less`. There should also be 17 chromosomes (if your file has fewer, try re-downloading the file using FireFox).
   ```bash
   # create the BWA index for fast read alignment:
   bwa index Scerevisiae.fasta

   # create the samtools index for fast search:
   samtools faidx Scerevisiae.fasta

   # create a dict file for GATK:
   samtools dict Scerevisiae.fasta >Scerevisiae.dict
   ```
   
5. Download the reads with SRR number `SRR10178655` from the SRA and convert to FASTQ format:
   ```bash
   # Fetch the .sra container file from the SRA FTP:
   wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR101/SRR10178655/SRR10178655.sra

   # Use SRA Toolkit to extract the reads in FASTQ format:
   fastq-dump --gzip --split-files --defline-seq '@$sn/$ri' --defline-qual '+' SRR10178655.sra
   ```

6. Run FastQC on the FASTQ files and examine the report (see `fastqc --help` for a complete list of options).
    - How many read pairs are in included in the FASTQ file?
      
      **ANSWER**: 2,335,497 pairs. FastQC (run separately on each end) will report the number of pairs in the "Total Sequences" field.
    - How long are the reads?
      
      **ANSWER**: 151 bp reported in "Sequence length" field. This number is the average length if the input sequences are of non-equal lengths.
    - What quality encoding are the reads in? What is the quality offset?
      
      **ANSWER**: The quality encoding is Sanger / Illumina 1.9. Looking at the slides on quality encoding, this means it is 33-offset.
    - For which metrics are there warnings?
      
      **ANSWER**: for both Read 1 and Read 2: "Per base sequence content" and "Per sequence GC content"; for Read 2 only: "Overrepresented sequences"
    - Are there any over-represented sequences in the file? If so, what is it?
      
      **ANSWER**: Yes. "Illumina Single End PCR Primer 1" which has sequence "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTATTAGCCGGTGTAGATCT". In a real analysis setting, we should add this sequence to our adapters file.

7. Next, run the Trimmomatic adapter trimmer on the FASTQ files in "PE" mode, using this [adapter file](https://www.dropbox.com/s/tpmhcz24jluq97s/adapters.fa). What fraction of the data were discarded? **NOTE**: Trimmomatic is in it's own subdirectory `sofware/Trimmomatic-0.39/trimmomatic-0.39.jar`.
    
    **ANSWER**: 6.45% = 100.0 - 93.55% both surviving; Large fractions of reads thrown out usually indicates that the sequencing library is contaminated with adapter dimers. This library looks pretty good.
    ```bash
    # assuming I'm working in my /Users/jessenbredeson/pfb2019/workshops/NGS directory:
    java -Xmx500m -jar ./software/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
        -summary SRR10178655.summary SRR10178655_1.fastq.gz SRR10178655_2.fastq.gz \
	SRR10178655_1_passed.fastq SRR10178655_1_failed.fastq \
	SRR10178655_2_passed.fastq SRR10178655_2_failed.fastq \
	MINLEN:100 ILLUMINACLIP:../pfb2019/workshops/NGS/data/adapters.fa:2:30:10:2:keepBothReads

    # I get the following statistics output to screen:
    Input Read Pairs: 2335497 Both Surviving: 2184756 (93.55%) Forward Only Surviving: 135634 (5.81%) Reverse Only Surviving: 326 (0.01%) Dropped: 14781 (0.63%)
    ```

8. Align the reads to the genome sequence using BWA-MEM. Convert the file to BAM format, sort the BAM file, and index it (see lecture notes for how).
   
   **ANSWER**:
   ```bash
    # be sure to include the ReadGroup (@RG) option!
    bwa mem -t 4 -R '@RG\tID:SRR10178655\tSM:SRR10178655\tLB:SRR10178655\tPL:ILLUMINA' \
       ./Scerevisiae.fasta SRR10178655_1_passed.fastq SRR10178655_2_passed.fastq \
       | samtools view -bS - >SRR10178655.bam

    # sort the BAM by coordinate using 1G of memory:
    samtools sort -m 1g -o SRR10178655.sorted.bam SRR10178655.bam

    # Be sure to index your BAM file!
    samtools index SRR10178655.sorted.bam
    ```

9. Now, Run GATK's `MarkDuplicates` tool on the BAM file to identify optical and PCR duplicate reads.
    
    **ANSWER**:
    ```bash
    gatk MarkDuplicates --java-options '-Xmx1G'  --VALIDATION_STRINGENCY SILENT \
        -I SRR10178655.sorted.bam -O SRR10178655.sorted.mdup.bam -M SRR10178655.sorted.mdup.metrics \
	-MAX_FILE_HANDLES 2000
    ```

10. Run `samtools stats` and `plot-bamstats` on the BAM file and examine the results.
    
    **ANSWER**: The commands to run:
    ```bash
    samtools stats --ref-seq Scerevisiae.fasta SRR10178655.sorted.mdup.bam >SRR10178655.stats
    plot-bamstats -s Scerevisiae.fasta >Scerevisiae.gc
    plot-bamstats -r Scerevisiae.gc -p SRR10178655 SRR10178655.stats
    open -a Safari.app SRR10178655.html
    ```
    - What is the mode insert size of the sequencing library?
        
	  **ANSWER**: Find the "Insert Size" plot (top-left plot). It should be ~252 bp.
	
    - What is the estimated read base-call error rate?
        
	  **ANSWER**: 0.77%, but this may be an over-estimate, as real polymorphisms count as mismatches
	
    - What fraction of your reads are duplicates?
        
	  **ANSWER**: 2.0%, high duplication percentages may be indicative of a low-complexity sequencing library, biased over-amplification of short fragments during PCR amplification, or the library was under-loaded onto the flowcell.
	
    - Within which base quality score range do the majority of mismatches occur? (*HINT*: see the "Mismatches per cycle" plot). Record the upper value of the range to use as the minimum base quality threshold for variant calling a later step.
        
	  **ANSWER**: From the "Mismatches per cycle" plot, you should observe that the yellowy-brown line representing `20 >= Q > 10` has an over-representation of base-call errors, so we want to set our minimum base-quality score to 20 to exclude them. To be really thorough, we should also check the alignments using `samtools tview` that these mismatches are in fact errors and not single-nucleotide variants.

11. Use `samtools view` to keep alignments with `PAIRED` and `PROPER_PAIR` flags *AND DO NOT* contain `UNMAP`, `MUNMAP`, `SECONDARY`, `QCFAIL`, `DUP`, or `SUPPLEMENTARY` flags; write the output as a BAM file. Index this new file with SAMtools (*HINT*: see `samtools flags` for help with flags). Be sure to index your filtered BAM file.
    
    **ANSWER**: We can create SAM/BAM flag-based filters by using `samtools flags` and combining the filters with commas:
    ```
    # Combine the binary array values for flags we DO want our reads to have:
    samtools flags PAIRED,PROPER_PAIR
    # returns: 0x3  3  PAIRED,PROPER_PAIR

    # Combine the binary array values for flags we DO NOT want our reads to have:
    samtools flags UNMAP,MUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY
    # returns: 0xf0c  3852  UNMAP,MUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY

    # use `samtools view` to output a new BAM file with our above filteres applied:
    samtools view -b -f 3 -F 3852 SRR10178655.sorted.mdup.bam >SRR10178655.sorted.mdup.filtered.bam
    ```
    
12. Using the base quality score determined in *Step 10* as the minimum base quality threshold, call SNPs using the GATK HaplotypeCaller. **NOTE**: using `--java-options '-Xmx1G'` below allows us to increate the maximum amount of memory allocated to GATK to run. If you ever encounter out-of-memory errors, increase this value (spaces not allowed between `-Xmx` and requested memory size).
    
    **ANSWER**:
    ```bash
    gatk HaplotypeCaller --java-options '-Xmx1G' \
        --minimum-mapping-quality 30 --min-base-quality-score 20 \
	--read-validation-stringency SILENT --reference Scerevisiae.fasta \
	--input SRR10178655.sorted.mdup.filter.bam --output SRR10178655.vcf
    ```

13. Use the `samtools depth` command to calculate the per-site depth of reads in the genome (see `samtools depth --help` for more info). The output file contains three columns: the chromosome name, position (1-based), and depth. For example:
    ```
    chrI	1	15
    chrI	2	15
    chrI	3	14
    chrI	4	16
    chrI	5	16
    ```
    
    **ANSWER**: To get the best results, also set the base quality threshold (`-q`) and mapping quality threshold (`-Q`) flags to the same as you used for variant calling:
    ```bash
    samtools depth -q 20 -Q 30 SRR10178655.sorted.mdup.filter.bam >SRR10178655.sorted.mdup.filter.depth
    ```
    
14. Write a script that computes a text histogram of depth with output similar to the following:
    ```
     1 |                                        
     2 |]                                       
     3 |]]                                      
     4 |]]]]]                                   
     5 |]]]]]]]]                                
     6 |]]]]]]]]]]]]]                           
     7 |]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]  
     8 |]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
     9 |]]]]]]]]]]]]]]]]]]]]]]]]]]]             
    10 |]]]]]]]]]]]]]]]                         
    11 |]]]]]]]]]                               
    12 |]]]]]                                   
    13 |]                                       
    14 |]                                       
    15 |                                        
    ```
    
    **ANSWER**: See [hist.py](hist.py); I chose min depth theshold of 20 and max depth threshold of 56.

15. Filter SNPs and Indels for variant loci within the center 95% of the depth distribution (use your distribution from above; estimating by eye is fine). To filter loci, use VCFtools:
    ```bash
    vcftools --recode --recode-INFO-all --stdout --max-missing 1 \
      --minDP <lower-threshold> --maxDP <upper-threshold> \
      --vcf <your.vcf> >your.filtered.vcf
    ```

16. Finally, how many of these SNPs and Indels intersect CDS features? (*HINT*: extract CDS features into a new GFF3 file and use `bedtools intersect` to do this to extract unique SNP loci).
    
    **ANSWER**: With my depth thresholds, I get 6192 variant loci intersecting CDS sequences:
    ```bash
    # Extract just the CDS features:
    grep -E '^#|CDS' Scerevisiae.gff3 >cds.gff3

    # perform the intersection, reporting each variant in the VCF uniquely (only once):
    bedtools intersect -nonamecheck -u -b cds.gff3 -a SRR10178655.filter.vcf | wc -l  # returns: 6192
    ```

**BONUS:** Try loading the genome FASTA, annotation GFF3, and filtered SNPs into [IGV](https://software.broadinstitute.org/software/igv)
   - Open IGV
   - Navigate through `Genomes` => `Load Genomes from File...` and select the genome FASTA file.
   - Navigate `File` => `Load from File...` and load a VCF/GFF3/BED file.


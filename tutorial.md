Tutorial
========

In this tutorial you will learn how to process the sequencing reads of a
single-cell prepared with the C1-CAGE protocol. It includes the
extraction of the UMIs and ribosomal RNAs with Tagdust2, the alignment
of the reads to a reference genome with bwa and the quantification of
single transcript molecules with umicount and the level1 script.

This tutorial assumes that all the software mentioned in the above
paragraph are already installed on your system and ready tu use from the
command-line. If this is not the case, please refer to the
'Prerequisite' section.

Background
----------

The CAGE [1] and CAGEscan [2] reads are output in compressed FASTQ
format by the MiSeq or HiSeq sequencer. If the multiplexing is done
entirely with Illumina's indexing system, then the reads are alread
demultiplexed after basecalling. 

In the following example we refer to the CAGE and CAGEscan reads as R1
and R2 (100\_S100\_L001\_R1\_001.fastq.gz and
100\_S100\_L001\_R2\_001.fastq.gz) respectively

Create output folders
---------------------

Start by creating all the folders where the output of each software will
be stored with mkdir.

    mkdir tagdust_r1 unzip_r2 extracted_reads cleaned_reads sai sampe genome_mapped properly_paired cagescan_pairs cagescan_frags level1

Extract UMIs
------------

Tag extraction is done with [TagDust2] (http://sourceforge.net/projects/tagdust/files/tagdust-2.13.tar.gz/download) v2.13. 
During this step, the unique molecular identifiers (UMI) are transferred from the sequence to
the read name, and the spacers are removed. Reads that do not
follow the architecture [3] are discarded.

    tagdust -t8 -o tagdust_r1/100_S100_L001_R1_001 -1 F:NNNNNNNN -2 S:TATAGGG -3 R:N ./test_data/100_S100_L001_R1_001.fastq.gz

    gunzip -c test_data/100_S100_L001_R2_001.fastq.gz > unzip_r2/100_S100_L001_R2_001.fq

After single-end mode extraction, the CAGEscan reads are then filtered
with the program [syncpairs](https://github.com/mmendez12/sync_paired_end_reads) to restore the pairing with the CAGE
reads. :

    syncpairs tagdust_r1/100_S100_L001_R1_001.fq unzip_r2/100_S100_L001_R2_001.fq extracted_reads/100_S100_L001_R1_001.fq extracted_reads/100_S100_L001_R2_001.fq

The reads are then filtered against the sequences of ribosomal genes
and synthetic RNA spikes using TagDust2 v2.13 in paired-end mode. 

This step requires to add the reads architecture in an external file. Here we call it SimpleArchitecture.txt and it contains an empty architecure. 
You can generate this file with the following command:

    echo 'tagdust -1 R:N' > 'SimpleArchitecture.txt'

Then you can filter the reads with:

    tagdust -arch SimpleArchitecture.txt -ref ercc_and_hg38_rRNA.fa -o cleaned_reads/100_S100_L001 extracted_reads/100_S100_L001_R1_001.fq extracted_reads/100_S100_L001_R2_001.fq

The file `ercc_and_hg38_rRNA.fa` contains the sequence of the spikes and the
human ribosomal RNA locus (GenBank ID U13369.1).  The sequence of the External
RNA Control Consortium (ERCC) spikes is in the public domain (see NIST's
[Certificate of Analysis](https://www-s.nist.gov/srmors/view_cert.cfm?srm=2374)
for [SRM 2374](https://www-s.nist.gov/srmors/view_detail.cfm?srm=2374)) and can
be [downloaded from NIST](https://www-s.nist.gov/srmors/view_datafiles.cfm?srm=2374).


Align paired-end reads
----------------------

For each sample separately, the reads are aligned paired-end on the
selected genome with [BWA](https://github.com/lh3/bwa) sampe  using standard parameters, except
"maximum\_insert\_size" that is set to 2,000,000. :

    bwa aln ./hg19_female.fa cleaned_reads/100_S100_L001_READ1.fq > sai/100_S100_L001_R1_001.sai
    bwa aln ./hg19_female.fa cleaned_reads/100_S100_L001_READ2.fq > sai/100_S100_L001_R2_001.sai

    bwa sampe -a 2000000 -c 0.00001 ./hg19_female.fa sai/100_S100_L001_R1_001.sai sai/100_S100_L001_R2_001.sai cleaned_reads/100_S100_L001_READ1.fq cleaned_reads/100_S100_L001_READ2.fq > sampe/100_S100_L001.sam

The alignments are converted to BAM format and sorted by coordinates
using [samtools](https://github.com/samtools/samtools/releases/latest). The result is the "genome mapped" reads:

    samtools view -uSo - sampe/100_S100_L001.sam | samtools sort - genome_mapped/100_S100_L001

CAGEscan fragments
------------------

Using samtools, the "genome mapped" reads are filtered by removing
non-properly paired reads [4] and non-primary alignments, and then sorted by
name. The result is the "properly paired" reads:

    samtools view -f 0x0002 -F 0x0100 -uo - genome_mapped/100_S100_L001.bam | samtools sort -n - properly_paired/100_S100_L001

The "properly paired" reads are then converted to BED12 format with the
program [pairedBamToBed12] (https://github.com/Population-Transcriptomics/pairedBamToBed12). These are the "CAGEscan pairs". :

    pairedBamToBed12 -i properly_paired/100_S100_L001.bam > cagescan_pairs/100_S100_L001.bed

The CAGEscan pairs sharing the same TSS and UMI are then aggregated into
"CAGEscan fragments", with the "umicountFP" script of [umicount] (https://github.com/mmendez12/umicount/).
The CAGEscan fragments represent single transcript molecules. :

    umicountFP -f cagescan_pairs/100_S100_L001.bed > cagescan_frags/100_S100_L001.bed

level1
------

An expression table of transcript counts is prepared with the "properly
paired" reads using the "level1" command of the 
[PromoterPipeline] (http://genome.gsc.riken.jp/plessy-20150516/PromoterPipeline_20150516.tar.gz)
version 2015.02.12 or greater (for UMI support), after filtering out CAGEscan reads with
the flag `0x40`. Reads that have both the same TSS and UMI count for one
transcript.

    python ./PromoterPipeline_20150516/level1.py -o level1/mylevel1file.l1.osc.gz -f 0x0042 -F 0x0104 --fingerprint genome_mapped/100_S100_L001.bam

This data has single nucleotide resolution, and can be converted
in to a promoter expression table with the "level2" command of the
PromoterPipeline, by distance based clustering with a default distance
of 20 between TSS.  Note that the input must be sorted, otherwise, the
PromoterPipeline scripts will produce incorrect data.

References/Notes
----------------

[1] The CAGE read is the 5′ sequence of the cDNA. On current Illumina
sequencers, it corresponds to “Read 1”.

[2] The CAGEscan read is the reverse-complement of the 3′ end of the
fragmented cDNA. On current Illumina sequencers, it corresponds to “Read
2”.


[3] The architecture is defined in TagDust2 in semantic blocks
that the program uses to build hidden markov models. In the following
example we use the following:

    tagdust -1 F:NNNNNNNN -2 S:TATAGGG -3 R:N


[4] The removal of "improper" pairs discards potentially valuable
information, if some transcripts were trans-splicing products, or if the
cells had chromosomal translocations that are not reflected in the
reference genome.





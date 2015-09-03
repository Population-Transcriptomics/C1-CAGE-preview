Prerequisites
=============

Here is a list of software to install and where to download them.

In the tutorial section we assume that the software are available from the command-line so you have 
to make you set you PATH correctly.

BWA
---

    git clone https://github.com/lh3/bwa.git
    cd bwa; make

After installing make sure to generate an index of the reference genome you need to use

    bwa index my_genome.fastq

tagdust
-------

    wget http://sourceforge.net/projects/tagdust/files/tagdust-2.13.tar.gz
    tar xzf tagdust-2.13.tar.gz
    cd tagdust-2.13
    ./configure
    make install
    
samtools
--------

    wget http://downloads.sourceforge.net/project/samtools/samtools/1.2/samtools-1.2.tar.bz2
    tar xjf samtools-1.2.tar.bz2
    cd samtools-1.2
    make install

pairedBamToBed12
----------------

    git clone https://github.com/Population-Transcriptomics/pairedBamToBed12.git
    cd pairedBamToBed12
    make

syncpairs
---------

    git clone https://github.com/mmendez12/sync_paired_end_reads.git
    cd sync_paired_end_reads
    python setup.py install

umicount
--------

    git clone https://github.com/mmendez12/umicount.git
    cd umicount
    python setup.py install


PromoterPipeline
----------------

    wget http://genome.gsc.riken.jp/plessy-20150516/PromoterPipeline_20150516.tar.gz
    tar xzf PromoterPipeline_20150516.tar.gz

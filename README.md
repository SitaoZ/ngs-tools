# Overview

The ngs-tools is python package for bioinformatics. Package format is a bioinformatics file format. Package util is a common tools software. Package PopGen is a population genetics software

1.1 **Install** <br>
    ```bash
    tar -zxvf Python-3.6.1.tgz
    cd Python-3.6.1
    ./configure
    make && make install
    git clone https://github.com/zhusitao1990/ngs-tools.git
    ```
1.2 **Fasta usage** <br>
   ```Python
    >>> from fasta import Fasta
    >>> fasta_obj = Fasta('/Users/zhusitao/data/chr1.fa')
    >>> len(fasta_obj)                      # get fasta object length
    >>> a_fasta = fasta_obj['a']            # index
    >>> fasta_dict = fasta_obj.readFasta()  # put fasta to a dict contain id and seq
    >>> fasta_ids = fasta_obj.fasta_key()   # return id list
    >>> fasta_seq = fasta_obj.fasta_sequence() # return seq list
    >>> fasta_obj.gc_rate("gc_rate.txt")       # return a gc_rate file
    >>> fasta_max_seq = fasta_obj.extract_item('max',"max_fasta.fa") # exact a max length fasta
    >>> fasta_min_seq = fasta_obj.extract_item('min',"min_fasta.fa") # exact a min length fasta
    ```

1.3 **Fastq usage** <br>
    ```Python
    >>> from fastq import Fastq
    >>> fastq_obj = Fastq('/Users/zhusitao/data/kio.fastq')
    >>> quality_system = fastq_obj.qualitySystem()   # return fastq quality system
    >>> fastq_obj.to_fasta("fasta_path")             # transform fastq to fasta
    >>> fastq_dict = fastq_obj.fastq_to_dict()       # return a dict dict[key] = fastq_record
    >>> index_seq = fastq_obj.indexSequence()        # return index sequence if exist
    >>> pairOrSingel = fastq_obj.pairOrSingel()      # return pair or singel
    ```

1.4 **GFF usage**
    ```Python
    >>> from gff import GFF
    >>> gff_obj = GFF("gff.path")
    >>> gff_dict = gff_obj.all_gff()      # return a dict[rownum] = gff_record
    >>> gene_number = gff_obj.geneCount() # return gene number
    >>> mRNA_number = gff_obj.mrnaCount() # return mRNA number
    >>> exon_number = gff_obj.exonCount() # return exon number
    >>> gff_obj.mRNA_fasta("mRNA_fa path")# exact mRNA fasta
    >>> gff_obj.cds_fasta("cds_fa path")  # cds fasta path
    ```

1.5 **fish usage**
    ```bash
    python3 fish.py -h

    usage: fish.py [-h] [-bf BAIT_FORMAT] [-ff FISH_FORMAT] [-bc BAIT_COLUMN]
               [-fc FISH_COLUMN] [-E] [--version]
               bait fish

    The program use bait to get wanted fish!

    positional arguments:
      bait                  bait file [table,fasta,gff,fq]
      fish                  fish file [table,fasta,gff,fq]

    optional arguments:
      -h, --help            show this help message and exit
      -bf BAIT_FORMAT, --bait_format BAIT_FORMAT
                        bait file formate [table,fasta,gff,fq]
      -ff FISH_FORMAT, --fish_format FISH_FORMAT
                        fish file formate [table,fasta,gff,fq]
      -bc BAIT_COLUMN, --bait_column BAIT_COLUMN
                        bait column for extract, default 1st column [-bc 1 or
                        --bait_column 1]
      -fc FISH_COLUMN, --fish_column FISH_COLUMN
                        fish column for extract, default 1st column [-fc 1 or
                        --fish_column 1]
      -E, --Except          if -e means contain in fish pool,or not
      --version             show program's version number and exit
    ```
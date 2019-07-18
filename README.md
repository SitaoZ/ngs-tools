# ngs-tools

Python program for bioinformatics. Package format is a bioinformatics file format. Package util is a common tools software. Package PopGen is a population genetics software

1. **Install** <br>
   git clone https://github.com/zhusitao1990/ngs-tools.git <br>

2. **Fasta example** <br>
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

3. **Fastq example** <br>
    ```Python
    >>> from fastq import Fastq
    >>> fastq_obj = Fastq('/Users/zhusitao/data/kio.fastq')
    >>> quality_system = fastq_obj.qualitySystem()   # return fastq quality system
    >>> fastq_obj.to_fasta("fasta_path")             # transform fastq to fasta
    >>> fastq_dict = fastq_obj.fastq_to_dict()       # return a dict dict[key] = fastq_record
    >>> index_seq = fastq_obj.indexSequence()        # return index sequence if exist
    >>> pairOrSingel = fastq_obj.pairOrSingel()      # return pair or singel
    ```


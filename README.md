# ngs-tools

Python program for bioinformatics. Package format is a bioinformatics file format. Package util is a common tools software. Package PopGen is a population genetics software

1. Install
   git clone https://github.com/zhusitao1990/ngs-tools.git

2. Example
   from fasta import Fasta
   fasta_obj = Fasta('/Users/zhusitao/data/chr1.fa')
   len(fasta_obj)
>>> fasta_dict = fasta_obj.readFasta()
>>> fasta_ids = fasta_obj.fasta_key()
>>> fasta_seq = fasta_obj.fasta_sequence()
>>> fasta_obj.gc_rate("gc_rate.txt")
>>> fasta_max_seq = fasta_obj.extract_item('max',"max_fasta.fa")
>>> fasta_min_seq = fasta_obj.extract_item('min',"min_fasta.fa")
>>> a_fasta = fasta_obj['a']

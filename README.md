# ngs-tools

Python program for bioinformatics. Package format is a bioinformatics file format. Package util is a common tools software. Package PopGen is a population genetics software

1. Install
   git clone https://github.com/zhusitao1990/ngs-tools.git

2. Example
   \>>> from fasta import Fasta <br>
   fasta_obj = Fasta('/Users/zhusitao/data/chr1.fa') <br>
   len(fasta_obj)                     # get fasta object length <br>
   fasta_dict = fasta_obj.readFasta() # put fasta to a dict contain id and seq <br>
   fasta_ids = fasta_obj.fasta_key()  # <br>
   fasta_seq = fasta_obj.fasta_sequence() <br>
   fasta_obj.gc_rate("gc_rate.txt") <br>
   fasta_max_seq = fasta_obj.extract_item('max',"max_fasta.fa") <br>
   fasta_min_seq = fasta_obj.extract_item('min',"min_fasta.fa") <br>
   a_fasta = fasta_obj['a'] <br>

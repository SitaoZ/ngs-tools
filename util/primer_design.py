# coding:utf-8

import re
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def gc(sequence):
    """ GC content accumulate
    # Arguments:
        sequence: a string of ATGC
    # Output
        a gc rate of sequence
    """
    gc_contend = (sequence.count('G') + sequence.count('C')) / len(sequence)
    return gc_contend


def primer_parser(seq, min_length=20, max_length=28):
    """used for primer parser
    # Arguements:
        seq: a string of ATGC
        min_length: min length for a primer
        max_length: max length for a primer
    # Output
        return a specific primer with direction 5' to 3'
    """
    try:
        if len(seq) < 20:
            raise Exception("sequence length less then 20 ")
        else:
            for i in range(min_length, max_length + 1):
                gc_content = gc(seq[:i])
                if 0.4 <= gc_content <= 0.6:
                    return seq[:i]
                else:
                    continue
            if gc(seq[:25]) < 0.4:
                return seq[:25]+"|low_abnormal"
            elif gc(seq[:25]) > 0.6:
                return seq[:22]+"|high_abnormal"

    except ValueError as error:
        print(error, 'sequence should be in 20-28')


def primer(seq, initiation, termination, direction, length=50):
    """used for primer design of sequence
    # Arguments
        seq : a Seq object from Bio.Seq with ATGC and a
            transcript contains 5'utr cds and 3'utr
        initiation : cds start with a ATG codon
        termination : cds end with a TAG/TAA/TGA codon
        direction : an integer 5 or 3
        length : primer length
    # Output
        return primer forward and primer reverse
        primer should be  5' to 3' direction
    """
    utr5 = seq[:initiation]
    cds = seq[initiation:termination]
    utr3 = seq[termination:]
    if direction == 3:
        reverse_pre = cds[-length:].reverse_complement()
        # 20-28 bp, GC 40% and 60%，contain TGA/TAG/TAA
        reverse = primer_parser(reverse_pre)
        return reverse
    elif direction == 5:
        reverse_pre = utr5[-length:].reverse_complement()
        # 20-28 bp, GC 40% and 60%，not contain ATG
        reverse = primer_parser(reverse_pre)
        return reverse


pattern_c = re.compile(r'CDS=(\d+?)-(\d+?)$')
pattern_t = re.compile(r'transcript:(\S+?)(\s)')

transcript = 'Ath_transcripts.fa'  # all transcript fasta from gffread in ensembl 39
abnormal_trans = []
five_primer = []
three_primer = []
a = 0
for record in SeqIO.parse(transcript, 'fasta'):
    chrom = record.id
    seq = record.seq
    desc = record.description
    matched_c = pattern_c.search(desc)
    matched_t = pattern_t.search(desc)
    a += 1
    if matched_c and matched_t:
        start, end = map(int, [matched_c.group(1), matched_c.group(2)])
        start = start - 1  # python list begin with 0
        transcript_id = matched_t.group(1)
        utr5, cds, utr3 = seq[:start], seq[start:end], seq[end:]
        if cds.startswith("ATG"):
            if 600 > len(utr5) > 50:
                # 5‘UTR 80% 343; 90% 480
                reverse = primer(seq, start, end, direction=5, length=50)
                rec = SeqRecord(
                    reverse,
                    id=transcript_id,
                    description="Arabidopsis thaliana 5utr_length=%d" % (len(utr5)))
                five_primer.append(rec)
            if 600 > len(utr3) > 50:
                # 3‘UTR 80% 408; 90% 554
                forward = primer(seq, start, end, direction=3, length=20)
                rec = SeqRecord(
                    forward,
                    id=transcript_id,
                    description="Arabidopsis thaliana 3utr_length=%d" % (len(utr3)))
                three_primer.append(rec)
        else:
            abnormal_trans.append(transcript_id)  # abnormal transcript
        if not a:
            break
SeqIO.write(five_primer, 'Ath_5utr_primer.fa', 'fasta')
SeqIO.write(three_primer, 'Ath_3utr_primer.fa', 'fasta')

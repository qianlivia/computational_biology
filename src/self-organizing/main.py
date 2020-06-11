from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped
from tools.self_organizing_tree import SOTA
import time
import numpy as np


def main():
    file_name = "data/coding.fa"
    # file_name = "data/cons_noncode.fa"
    alignment = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
    for seq_record in SeqIO.parse(file_name, "fasta"):
        alignment.extend([seq_record])

    par = SOTA(alignment[10:15])
    start = time.time()
    par.train()
    end = time.time()
    print("Self-organizing tree network ran in {} seconds.".format(end - start))
    par.draw_tree()

if __name__ == '__main__':
    main()
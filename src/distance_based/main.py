from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped
from distance_calculator import Distance_Calculator

file_name = "data/coding.fa"
alignment = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
for seq_record in SeqIO.parse(file_name, "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    alignment.extend([seq_record])

# Print the alignment
print(alignment)

dc = Distance_Calculator()
dm = dc.create_distance_matrix(alignment)
dm.data.to_csv("animals.csv")
dc.build_tree(dm)
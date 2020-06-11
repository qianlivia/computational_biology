from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped
import time

# Prettify labels
def get_label(leaf):
    if leaf.name.startswith("Inner"):
        return ""
    return leaf.name.replace("_", " ")

# Read the sequences and align

aln = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
for seq_record in SeqIO.parse("data/coding.fa", "fasta"):
# for seq_record in SeqIO.parse("data/cons_noncode.fa", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    aln.extend([seq_record])

# Print the alignment
print(aln)

# Calculate the distance matrix
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)

# Print the distance Matrix
print('\nDistance Matrix\n===================')
print(dm)

# Construct the phylogenetic tree using UPGMA algorithm
constructor = DistanceTreeConstructor()

print("Neighbor joining")
start = time.time()
tree = constructor.nj(dm)
end = time.time()
print("Neigbor joining ran in {} seconds.".format(end - start))

# Draw the phylogenetic tree
Phylo.draw(tree, label_func=get_label)

# Print the phylogenetic tree in the terminal
print('\nPhylogenetic Tree\n===================')
Phylo.draw_ascii(tree)

print("----------------------")

"""
print("UPGMA")

start = time.time()
tree = constructor.upgma(dm)
end = time.time()
print("UPGMA ran in {} seconds.".format(end - start))

# Draw the phylogenetic tree
Phylo.draw(tree, label_func=get_label)

# Print the phylogenetic tree in the terminal
print('\nPhylogenetic Tree\n===================')
Phylo.draw_ascii(tree)

print("----------------------")
"""
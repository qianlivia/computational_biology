from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Alphabet import IUPAC, Gapped
from Bio import Phylo
import time

def get_label(clade):
    if clade.name.startswith("Inner"):
        return ""
    return clade.name.replace("_", " ").split("/")[0]

def main():
    file_name = "data/coding.fa"
    # file_name = "data/cons_noncode.fa"

    alignment = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
    for seq_record in SeqIO.parse(file_name, "fasta"):
        alignment.extend([seq_record])

    print("Number of characters in alignment:", len(alignment[0]))

    ####################
    # Neighbor joining #
    ####################
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    start = time.time()
    tree = constructor.nj(dm)
    end = time.time()
    print("Neighbor joining ran in {} seconds.".format(end - start))
    Phylo.draw(tree, label_func=get_label)

    #########
    # UPGMA #
    #########

    start = time.time()
    tree = constructor.upgma(dm)
    end = time.time()
    print("UPGMA ran in {} seconds.".format(end - start))
    Phylo.draw(tree, label_func=get_label)

if __name__ == '__main__':
    main()
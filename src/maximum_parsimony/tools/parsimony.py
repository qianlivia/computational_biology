import numpy as np
from Bio.Phylo import BaseTree
from Bio.Phylo import draw
from Bio.Align import MultipleSeqAlignment

class Parsimony():

    def __init__(self, alignment: MultipleSeqAlignment):
        # self.tree = BaseTree.Tree(None, rooted=False)
        self.alignment = alignment.copy()
        self.size = len(alignment._records)
        self.alignment_length = alignment.get_alignment_length()
        self.inner_nodes_parents = np.zeros((self.size - 1))
        self.leaves_parents = np.zeros((self.size))

    def calc_parsimony(self):
        """
        Calculate score for every site (self.alignment_length).
        """
        pass

    def parsimony_for_site(self, tree: BaseTree, u: int):
        """
        Calculate parsimony score for given site. From bottom to the top.
        https://www.cs.helsinki.fi/bioinformatiikka/mbi/courses/07-08/itb/slides/itb0708_slides_158-191.pdf
        https://www.bio.fsu.edu/~stevet/BSC5936/Swofford.F2003.2.pdf
        https://tel.archives-ouvertes.fr/tel-01479049/document

        u: given site
        """
        pass
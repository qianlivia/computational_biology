import numpy as np
from math import ceil
from Bio.Phylo import BaseTree, draw
from Bio.Phylo.BaseTree import Clade
from Bio.Align import MultipleSeqAlignment
from enum import Enum

class Base(Enum):
    """
    Enum for different bases and their index values.
    """
    A = 0
    C = 1
    G = 2
    T = 3


class SOTA():
    """
    Self-Organizing Tree Algorithm. The algorithm presented here is based both on the Kohonen self-organizing maps
    and on the growing cell structures algorithm of Fritzke.
    Based on a publication by J. Dopazo (1997), doi: 10.1007/PL00006139.
    """

    def __init__(self, alignment: MultipleSeqAlignment, eta: double = 0.2):
        """
        Constructor.
        
        alignment: MultipleSeqAlignment containing the alignment
        """

        self.size = len(alignment) # Number of alignments
        self.length = alignment.get_alignment_length() # Length of each sequence

        # If the number of sequences isn't enough
        if self.size <= 0:
            raise ValueError("There aren't enough taxa.")

        # Initialize tree with 1 node. To be continuously updated.
        inner_clade = ParsimonyClade(None, None, sets=[], score=0)
        first_cell = ParsimonyClade(None, None, sets=[], score=0)
        inner_clade.clades.append(first_cell)
        self.tree = inner_clade

        self.cells = [first_cell] # Cells in order of creation

        # Sequences to classify. Dimension: number of taxa, number of different bases + 1, length of sequences
        self.S = np.zeros((self.size, len(Base) + 1, self.length))
        self.names = [] # Names of animal species belonging to the sequences.
        
        # Weights belonging to cells. Dimension: number of initial cells, number of different bases + 1, length of sequences
        self.C = np.random.uniform(0, 1, (1, len(Base) + 1, self.length))

        # Mapping telling which cell each sequence belongs to (element in position i should be assigned to cell mapping[i]).
        # Subject to continuous update.
        self.mapping = [0] * self.size

        # For every alignment, do the corresponding coding.
        for i in range(self.size):
            sequence = alignment[i]
            self.names.append(sequence.id) # Store name

            for j in range(self.length):
                char = sequence.seq[j]
                try:
                    base = Base[char]
                    self.S[i][base.value][j] = 1
                except KeyError:
                    self.S[i][len(Base) + 1][j] = 1

    def train(self, epochs: int = 20):
        """
        The algorithm.

        epochs: number of epochs for adaptation
        """

        # If there is only one taxon
        if self.size == 1:
            cell = self.tree.clades[0]
            self.tree = self.create_tree(self.tree)
            return

        # Return if there are only two taxa
        if self.size == 2:
            self.tree = self.create_tree(self.tree)
            return

        # Grow until the tree is complete (until it has number of leaves sequal to self.size).
        for i in range(self.size):

            # Step 1: Adaptation for a number of epochs.
            iters = 0
            while iters < epochs:
                for i in range(self.size):
                    # Fetch corresponding sequence.
                    s = self.S[i] # 2D matrix

                    # Min search: find the closest cell
                    min_sim_ind = self._find_most_similar_cell(s)
                    self.mapping[i] = min_sim_ind
                    
                    # TODO Weight update
                    # Find the index of neighbors

                # To the next iteration
                iters += 1

            # TODO: Step 2: Growing the network.

    def _find_most_similar_cell(self, s):
        """
        Find the most similar cell.

        s: sequence (matrix)
        returns: index of most similar cell
        """
        min_sim = np.inf
        min_sim_ind = -1
        for j in range(len(self.cells)):
            sim = self._similarity(s, self.C[j])
            if sim <= min_sim:
                min_sim = sim
                min_sim_ind = j
        return min_sim_ind

    def _similarity(self, s, c):
        # TODO: implement
        return 0

    def draw_tree(self, show_branch_labels: bool = False):
        """
        Draws tree.

        show_branch_labels: show branch lengths (optional); default is False
        """
        if self.tree is None:
            print("Please first build the tree.")
        else:
            if show_branch_labels:
                draw(self.tree, label_func=self.get_label, branch_labels=(lambda clade: clade.branch_length))
            else:
                draw(self.tree, label_func=self.get_label)
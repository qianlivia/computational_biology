import numpy as np
from Bio.Phylo import BaseTree, draw
from Bio.Phylo.BaseTree import Clade
from Bio.Align import MultipleSeqAlignment
from enum import Enum
from collections import Counter

class Base(Enum):
    """
    Enum for different bases and their index values.
    """

    A = 0
    C = 1
    G = 2
    T = 3


class SOTClade(Clade):
    """
    Class extending Bio.Phylo.BaseTree.Clade. Stores base sets and parsimony scores.
    """
    
    def __init__(
        self,
        branch_length=None,
        name=None,
        clades=None,
        confidence=None,
        color=None,
        height=1,
        width=1
    ):
        """
        Constructor.

        height: height of weight matrix
        width: width of weight matrix
        """
        
        super(SOTClade, self).__init__(branch_length, name, clades, confidence, color, width)
        
        # weight: weight belonging to the current cell
        self.weight = np.random.uniform(0, 1, (height, width))

    def copy(self):
        """
        Copy current instance.

        returns: new instance
        """

        obj = type(self).__new__(self.__class__)
        obj.__dict__.update(self.__dict__)
        return obj


class SOTA():
    """
    Self-Organizing Tree Algorithm. The algorithm presented here is based both on the Kohonen self-organizing maps
    and on the growing cell structures algorithm of Fritzke.
    Based on a publication by J. Dopazo (1997), doi: 10.1007/PL00006139.
    """

    def __init__(self, alignment: MultipleSeqAlignment):
        """
        Constructor.
        
        alignment: MultipleSeqAlignment containing the alignment
        """

        self.size = len(alignment) # Number of alignments
        self.length = alignment.get_alignment_length() # Length of each sequence
        self.nr_of_bases = len(Base) + 1
        self.alignment = alignment

        # If the number of sequences isn't enough
        if self.size <= 0:
            raise ValueError("There aren't enough taxa.")

        # If there is only one taxon
        if self.size == 1:
            inner_clade = Clade(None, None)
            first_cell = Clade(None, alignment[0].id)
            inner_clade.clades.append(first_cell)
            self.tree = self.create_tree(inner_clade)

        # If there are only 2 taxa
        elif self.size == 2:
            inner_clade = Clade(None, None)
            first_cell = Clade(None, alignment[0].id)
            second_cell = Clade(None, alignment[1].id)
            inner_clade.clades.append(first_cell)
            inner_clade.clades.append(second_cell)
            self.tree = self.create_tree(inner_clade)

        # In any other case
        else:
            self.reset()

            # Sequences to classify. Dimension: number of taxa, number of different bases + 1, length of sequences
            self.S = np.zeros((self.size, self.nr_of_bases, self.length))
            self.names = [] # Names of animal species belonging to the sequences.
            
            # For every alignment, do the corresponding coding.
            for i in range(self.size):
                sequence = self.alignment[i]
                self.names.append(sequence.id) # Store name

                for j in range(self.length):
                    char = sequence.seq[j]
                    try:
                        base = Base[char]
                        self.S[i][base.value][j] = 1
                    except KeyError:
                        self.S[i][self.nr_of_bases - 1][j] = 1

    def reset(self):
        """
        Reset variables.
        """

        # Initialize tree with 2 nodes. To be continuously updated.
        inner_clade = Clade(None, None)
        first_cell = SOTClade(None, None, height=self.nr_of_bases, width=self.length)
        second_cell = SOTClade(None, None, height=self.nr_of_bases, width=self.length)
        inner_clade.clades.append(first_cell)
        inner_clade.clades.append(second_cell)
        self.tree = inner_clade

        # Cells in order of creation (deletion included).
        # Dimension of weights belonging to cells: number of different bases + 1, length of sequences
        self.C = [first_cell, second_cell]

        self.parents = {} # Dictionary for storing the parents of nodes and cells
        self.parents[first_cell] = inner_clade # Store the first parent-child relationship
        self.parents[second_cell] = inner_clade

        # Mapping telling which cell each sequence belongs to (element in position i should be assigned to cell mapping[i]).
        # Subject to continuous update.
        self.mapping = {}

    def train(self, epochs: int = 20, eta: float = 0.1, alpha: list = [0.5, 0.5, 0.5], b = 0):
        """
        The algorithm.

        epochs: number of epochs for adaptation
        eta: learning rate
        """

        self.epochs = epochs
        self.eta = eta
        self.alpha = alpha
        self.b = b

        # There's more than one solution if there are more than 2 taxa.
        if self.size > 2:
            
            # Reset variables.
            self.reset()

            # Grow until the tree is complete (until it has number of leaves sequal to self.size).
            for i in range(self.size):

                # Step 1: Adaptation for a number of epochs.
                iters = 0
                while iters < epochs:
                    for i in range(self.size):
                        # Fetch corresponding sequence.
                        s = self.S[i] # 2D matrix - weight belonging to sequence i

                        # Min search: find the closest cell
                        cell = self._find_most_similar_cell(s)
                        self.mapping[i] = cell
                        
                        # Weight update
                        self._weight_update(s, cell.weight, self.alpha[0])

                        # Find the index of neighbor(s): a cell is considered a neighbor if they have a mutual parent. If the current
                        # cell (leaf) doesn not have a sister that's also a cell, only itself gets updated.
                        parent = self.parents[cell]

                        if parent.clades[0] == cell:
                            sister = parent.clades[1]
                        else:
                            sister = parent.clades[0]

                        # If sister is a cell, update (along with the mutual parent)
                        if sister.clades:
                            self._weight_update(s, parent.weight, self.alpha[1])
                            self._weight_update(s, sister.weight, self.alpha[2])

                    # To the next iteration
                    iters += 1

                # Step 2: Growing the network.
                # Find the cell that got associated with the highest number of sequences.
                cell, _ = Counter(self.mapping.values()).most_common(1)[0]

                # Create children.
                first_child = SOTClade(None, None)
                second_child = SOTClade(None, None)
                first_child.weight = cell.weight
                first_child.weight = cell.weight

                # Create relations.
                self.parents[first_child] = cell
                self.parents[second_child] = cell
                cell.clades.append(first_child)
                cell.clades.append(second_child)

                # Remove the current node from the list of cells.
                self.C.remove(cell)

                # Add new elements to the list.
                self.C.append(first_child)
                self.C.append(second_child)

            # Finalize tree by assigning labels.

            for i in range(self.size):
                cell = self.mapping[i]
                cell.name = self.names[i]

            self.tree = self.create_tree(self.tree)

    def _find_most_similar_cell(self, s: np.ndarray):
        """
        Find the most similar cell.

        s: sequence (matrix)
        returns: most similar cell
        """

        min_val = np.inf
        min_cell = None
        for cell in self.C:
            sim = self._similarity(s, cell.weight)
            if sim <= min_val:
                min_val = sim
                min_cell = cell
        return min_cell

    def _similarity(self, s: np.ndarray, c: np.ndarray):
        """
        Calculate distance between sequence s and cell c.

        s: sequence
        c: cell
        returns: distance
        """

        d = 0

        for l in range(self.length):

            temp = 0
            for r in range(self.nr_of_bases):
                temp += s[r][l] * c[r][l]

            d += 1 - temp

        d = d / self.length

        return d

    def _weight_update(self, s: np.ndarray, c: np.ndarray, alpha: float):
        """
        Update weights.
        
        s: array containing weight belonging to a sequence
        c: array containing weight belonging to a cell
        a: constant update rate depending on the role (winning cell, ancestor or sister cell)
        """

        t = self.size * self.epochs
        tau = self.epochs

        M = self.eta * self.nr_of_bases * self.length
        Eta = alpha * ((1 - t) / M) * (1 - self.b * tau)
        s += Eta * (s - c)

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
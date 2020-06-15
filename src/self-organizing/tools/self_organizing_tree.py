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


class SOTAClade(Clade):
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
        width=None,
        rows=1,
        columns=1
    ):
        """
        Constructor.

        height: height of weight matrix
        width: width of weight matrix
        """
        
        super(SOTAClade, self).__init__(branch_length, name, clades, confidence, color, width)
        
        # weight: weight belonging to the current cell
        self.weights = np.random.uniform(0, 1, (rows, columns))

    def set_weight(self, weights: np.ndarray):
        """
        Setter for weight matrix.

        weight: item to set
        """

        self.weights = weights

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

    Cycle: the series of operations performed until a cell generates two descendants.
    Presentation implies two steps:
        first, to find the best matching cell (winning cell) for each input sequence i and
        second, to update this cell and its neighborhood.
    """

    def __init__(self, alignment: MultipleSeqAlignment, seed: int = 0):
        """
        Constructor.
        
        alignment: MultipleSeqAlignment containing the alignment
        """

        self.seed = seed
        np.random.seed(seed)

        self.size = len(alignment) # Number of alignments
        self.length = alignment.get_alignment_length() # Length of each sequence
        self.nr_of_bases = len(Base) + 1
        self.alignment = alignment
        self.threshold = 10e-5
        self.E = 10e-5

        # If the number of sequences isn't enough
        if self.size <= 0:
            raise ValueError("There aren't enough taxa.")

        # If there is only one taxon
        if self.size == 1:
            inner_clade = SOTAClade(None, None)
            first_cell = SOTAClade(None, alignment[0].id)
            inner_clade.clades.append(first_cell)
            self.tree = self.create_tree(inner_clade)

        # If there are only 2 taxa
        elif self.size == 2:
            inner_clade = SOTAClade(None, None)
            first_cell = SOTAClade(None, alignment[0].id)
            second_cell = SOTAClade(None, alignment[1].id)
            inner_clade.clades.append(first_cell)
            inner_clade.clades.append(second_cell)
            self.tree = self.create_tree(inner_clade)

        # In any other case
        else:
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
                        self.S[i][base.value, j] = 1
                    except KeyError:
                        self.S[i][self.nr_of_bases - 1, j] = 1

    def reset(self):
        """
        Reset variables.
        """

        # Initialize tree with 1 node. To be continuously updated.
        inner_clade = SOTAClade(None, "Node0", rows=self.nr_of_bases, columns=self.length)
        self.tree = self.create_tree(inner_clade)
        self.nr_of_nodes = 1

        # Cells in order of creation (deletion included).
        # Dimension of weights belonging to cells: number of different bases + 1, length of sequences
        self.C = [inner_clade]

        # Dictionary for storing the parents of nodes and cells.
        self.parents = {}

        # Mapping telling which cell each sequence belongs to (element in position i should be assigned to cell mapping[i]).
        # Subject to continuous update.
        self.mapping = {}

        # Mapping telling which sequence belongs to a given cell.
        self.mapping_inv = {}
        self.mapping_inv[inner_clade] = []

    def train(self, epochs: int = 50, eta: float = 0.5, alpha: list = [0.5, 0.02, 0.01], b = 0):
        """
        The algorithm.

        epochs: number of epochs for adaptation
        eta: learning rate
        """

        self.epochs = epochs
        self.eta = eta
        self.alpha = alpha
        self.b = b

        # There's only one solution when the number of taxa is smaller than 3.
        if self.size <= 2:
            return
            
        # Reset variables.
        self.reset()

        # Grow until the tree is complete (until the resource threshold is reached).
        
        # Step 1: Adaptation.
        self._adaptation()
        _, res = self._find_highest_resource()

        it = 1
        while res > self.threshold:
            # Note: each iteration is a cycle.

            # Step 2: Growing the network.
            self._growing()

            # Step 3: Adaptation after growing.
            self._adaptation()
            
            # Find the cell with the highest resource value.
            _, res = self._find_highest_resource()

            it += 1

        # Finalize tree by assigning labels.
        for i in range(self.size):
            cell = self.mapping[i]
            cell.name = self.names[i]

    def _adaptation(self):
        """
        Adaptation process.
        """

        # Number of epochs dependent on error

        for iters in range(self.epochs):
            
            for i in range(self.size):
                # Note: each iteration of this loop is considered a presentation.

                # Fetch corresponding sequence.
                s = self.S[i] # 2D matrix - weight belonging to sequence i

                # Min search: find the closest cell
                cell = self._find_most_similar_cell(s)
                
                # Delete old mapping (if exists)
                try:
                    old_cell = self.mapping[i]
                    self.mapping_inv[old_cell].remove(i)
                except KeyError:
                    pass

                # Create new mapping between the sequence and the cell
                self.mapping[i] = cell
                self.mapping_inv[cell].append(i)
                
                # Weight update
                self._weight_update(s, cell.weights, self.alpha[0])

                # Find the index of neighbor(s): a cell is considered a neighbor if they have a mutual parent. If the current
                # cell (leaf) does not have a sister that's also a cell, only it itself gets updated.

                try:
                    parent = self.parents[cell]

                    if parent.clades[0] == cell:
                        sister = parent.clades[1]
                    else:
                        sister = parent.clades[0]

                    # If sister is a cell, update (along with the mutual parent)
                    if sister.clades:
                        self._weight_update(s, parent.weights, self.alpha[1])
                        self._weight_update(s, sister.weights, self.alpha[2])

                # Ignore if the cell is root.
                except KeyError:
                    pass

            if iters == 0:
                prev_err = self.get_error()
            else:
                new_err = self.get_error()
                relative_inc = np.abs((new_err - prev_err) / prev_err)
                if relative_inc < self.E or new_err < self.threshold:
                    break
                prev_err = new_err

    def _growing(self):
        """
        Growing of the network.
        """

        # Find the cell with the highest resource value.
        cell, res = self._find_highest_resource()

        # Create children.
        first_child = SOTAClade(None, "Node{}".format(self.nr_of_nodes))
        self.nr_of_nodes += 1
        second_child = SOTAClade(None, "Node{}".format(self.nr_of_nodes))
        self.nr_of_nodes += 1
        first_child.set_weight(cell.weights.copy())
        second_child.set_weight(cell.weights.copy())

        # Create relations.
        self.parents[first_child] = cell
        self.parents[second_child] = cell
        cell.clades.append(first_child)
        cell.clades.append(second_child)

        # Remove the current node from the list of cells.
        self.C.remove(cell)
        del self.mapping_inv[cell]

        # Add new elements to the list.
        self.C.append(first_child)
        self.mapping_inv[first_child] = []
        self.C.append(second_child)
        self.mapping_inv[second_child] = []

    def _find_highest_resource(self):
        """
        Find the cell with the highest resource value.
        """
        
        cell = self.C[0]
        res = self._get_resource(cell)

        for i in range(1, len(self.C)):
            c = self.C[i]
            r = self._get_resource(c)
            if r > res:
                res = r
                cell = c

        return cell, res

    def _find_most_similar_cell(self, s: np.ndarray):
        """
        Find the most similar cell.

        s: sequence (matrix)
        returns: most similar cell
        """

        min_val = np.inf
        min_cell = None
        for cell in self.C:
            sim = self._similarity(s, cell.weights)
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
        
        d = np.sum([1 - np.dot(s[:, l], c[:, l]) for l in range(self.length)])
        return d / self.length

    def _get_resource(self, cell: SOTAClade):
        """
        Get resurce value of cell.

        c: cell
        returns: resource value
        """

        count = 0
        summ = 0
        sequences = self.mapping_inv[cell]
        for seq in sequences:
            summ += self._similarity(self.S[seq], cell.weights)
            count += 1
        
        if count == 0:
            return 0
        return summ / count

    def _weight_update(self, s: np.ndarray, c: np.ndarray, alpha: float):
        """
        Update weights.
        
        s: array containing weight belonging to a sequence
        c: array containing weight belonging to a cell
        a: constant update rate depending on the role (winning cell, ancestor or sister cell)
        """

        Eta = alpha * self.eta
        c += Eta * (s - c)

    def get_error(self):
        """
        Get total error.

        returns: total error
        """

        summ = 0
        for seq, cell in self.mapping.items():
            summ += self._similarity(self.S[seq], cell.weights)

        return summ

    def create_tree(self, root: Clade):
        """
        Create tree with given root.

        root: root

        returns: tree
        """

        return BaseTree.Tree(root, rooted=True)

    def get_label(self, clade: SOTAClade):
        """
        Prettify tree labels.

        clade: current clade

        returns: the label belonging to the current clade
        """

        if clade.name is None:
            return ""
        if clade.name.startswith("Node"):
            return ""
        return clade.name.replace("_", " ").split("/")[0]

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
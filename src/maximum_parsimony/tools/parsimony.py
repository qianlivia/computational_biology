
from Bio.Phylo import BaseTree, draw
from Bio.Phylo.BaseTree import Clade
from Bio.Align import MultipleSeqAlignment
import numpy as np
from abc import abstractmethod
from collections import Counter

class ParsimonyClade(Clade):
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
        sets=None,
        score=0
    ):
        """
        Constructor.

        sets: list of sets containing bases (for each site in a sequence)
        score: parsimony score belonging to this clade
        """

        super(ParsimonyClade, self).__init__(branch_length, name, clades, confidence, color, width)
        self.sets = sets
        self.score = score


class Parsimony():
    """
    Class for methods using maximum parsimony.
    """

    def __init__(self, alignment: MultipleSeqAlignment):
        """
        Constructor.

        alignment: MultipleSeqAlignment containing the alignment
        """
        self.tree = None # Best tree
        self.trees = [] # Complete trees
        self.threshold = np.inf # The score of the currently best tree

        # Create the dictionary (and list) of taxa
        self.leaves = [] # Leaves as clades

        # For every alignment, create corresponding leaves and sets.
        for sequence in alignment:

            # Create a list of sets containing bases
            base_sets = []
            for char in sequence.seq:
                base_sets.append({char})

            # Create leaf
            leaf = ParsimonyClade(None, sequence.id, sets=base_sets, score=0)
            self.leaves.append(leaf)

        self.size = len(self.leaves) # Number of alignments
        self.length = alignment.get_alignment_length() # Length of each sequence

    @abstractmethod
    def run(self, print_best: bool = False):
        """
        Run the algorithm.

        print_best: print best tree (optional). Default is False.
        """
        pass

    def create_tree(self, root: ParsimonyClade):
        """
        Create tree with given root.

        root: root

        returns: tree
        """
        return BaseTree.Tree(root, rooted=False)

    def _create_inner_node(self, left: ParsimonyClade, right: ParsimonyClade, threshold: int, name: str = None):
        """
        Create new root with given left and right subtree, calculate the parsimony score and create the corresponding sets
        according to Fitch's algorithm. Return prematurely if the score exceeds the given threshold.

        left: root of left subtree
        right: root of right subtree
        threshold: the best result so far; the lower the better (optional). Default is infinity.
        name: name of new node (optional). Default is None.

        returns: whether the resulting parsimony score is under the given threshold; root of new tree (if first return value is true)
        """

        # Calculate the initial parsimony score by adding the scores of the left and the right node
        print(left, right)
        score = left.score + right.score
        # If this is already greater than or equal to the threshold, return
        if score >= threshold:
            return False, None

        # Create inner node and append left and right node
        root = ParsimonyClade(None, name)
        root.clades.append(left)
        root.clades.append(right)

        # Initialize variables
        sets = []
        
        # Count parsimony score for every site (self.length elements)
        for u in range(self.length):
            u_score, u_set = self._parsimony_for_site(left, right, u)
            score += u_score

            if score >= threshold:
                return False, None

            sets.append(u_set)

        # Set new score and sets
        root.score = score
        root.sets = sets

        return True, root

    def _parsimony_for_site(self, left: ParsimonyClade, right: ParsimonyClade, u: int):
        """
        Calculate the parsimony score of a given site using the sets and scores of children nodes
        (there's no need for access to the parent node).
        https://www.cs.helsinki.fi/bioinformatiikka/mbi/courses/07-08/itb/slides/itb0708_slides_158-191.pdf
        https://www.bio.fsu.edu/~stevet/BSC5936/Swofford.F2003.2.pdf
        https://tel.archives-ouvertes.fr/tel-01479049/document
        http://pages.stat.wisc.edu/~larget/Genetics629/Fall2009/outline2.pdf

        left: left child
        right: right child
        u: site (position in a sequence)
        """
        
        X = left.sets[u].intersection(right.sets[u]) # Intersection
        Y = left.sets[u].union(right.sets[u]) # Union
        score = 0 # Initial score
        
        # If not empty, X is the set belonging to the current node.
        if X:
            sets = X
        # Otherwise, the set for the node is the union of the children sets and one is added to the score.
        else:
            sets = Y
            score += 1

        return score, sets
    
    def get_label(self, clade: ParsimonyClade):
        """
        Prettify tree labels.

        clade: current clade

        returns: the label belonging to the current clade
        """
        if clade.name is None:
            return ""
        if clade.name.startswith("Inner"):
            return ""
        return clade.name.replace("_", " ")

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
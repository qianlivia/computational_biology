import numpy as np
from Bio.Phylo import BaseTree
from Bio.Phylo import draw
from Bio.Align import MultipleSeqAlignment
from collections import Counter
from tools.utils import remove_unprocessed, ParsimonyClade
import random

class Parsimony():
    """
    Class for methods using maximum parsimony.
    """

    def __init__(self, alignment: MultipleSeqAlignment, bnb: bool = True):
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
        self.bnb = bnb # Branch and bound. If false, use exhaustive search.
        random.shuffle(self.leaves)
        random.shuffle(self.leaves)

    def run(self, print_best = False):
        """
        Generate all possible trees and find the one with the lowest score.
        """
        if self.size <= 0:
            print("There aren't enough taxon.")
            return

        # If there is only one taxon
        if self.size == 1:
            inner_clade = ParsimonyClade(None, None, sets=[], score=0)
            inner_clade.clades.append(self.leaves[0])
            self.tree = self.create_tree(inner_clade)
            return
            
        # Initialize
        trees = [self.leaves[0]]

        # Depth first search
        self._DFS_trees(trees, 1)
        self.trees = [self.create_tree(tree) for tree in self.trees]
        self.tree = self.create_tree(self.tree)
        
        # Print trees
        for tree in self.trees:
            print(tree)
            print()

        print(len(self.trees))

        if print_best:
            print(self.tree)

    def create_tree(self, root: ParsimonyClade):
        """
        Create tree with given root.

        root: root

        returns: tree
        """
        return BaseTree.Tree(root, rooted=False)

    def _DFS_trees(self, queue: list, leaf_no: int):
        while queue:
            tree = queue.pop(0)

            if leaf_no < 8 - 1:
                # If this is not the last iteration
                new_trees = self._add_leaf(tree, self.leaves[leaf_no], threshold = self.threshold, final_iteration=False)
                self._DFS_trees(new_trees, leaf_no + 1)
            else:
                # If this is the last iteration
                new_trees = self._add_leaf(tree, self.leaves[leaf_no], threshold = self.threshold, final_iteration=True)
                # When using extensive search, the following are all the complete trees; when using BnB, these are the all the
                # temporarily best complete trees.
                self.trees.extend(new_trees)

    def _add_leaf(self, root: ParsimonyClade, leaf: ParsimonyClade, threshold: int = np.inf, final_iteration: bool = False):
        """
        Add leaf to a tree.

        root: root defining the given tree
        leaf: leaf to be added
        threshold: the best result so far (the lower the better)
        final_iteration: if this is the final iteration (when the last taxon needs to be added)

        returns: trees after the leaf has been added in every possible combination (represented by their roots)
        """

        new_trees = []
        min_tree = None

        ####################
        # Branch and bound #
        ####################
        if self.bnb:
            # There are three different kinds of tree formation.

            # No. 1.: create a new root; let the old root be the left subtree and the leaf the right one.
            under_threshold, inner_clade = self._create_inner_node(root, leaf, threshold = threshold)
            if under_threshold:
                if final_iteration:
                    threshold = inner_clade.score
                    min_tree = inner_clade
                new_trees.append(inner_clade)

            # If the root is a leaf, return current state (it doesn't have subtrees).
            if not root.clades:
                if final_iteration:
                    self.tree = min_tree
                    self.threshold = threshold
                return new_trees

            left = root.clades[0] # Left subtree
            right = root.clades[1] # right subtree

            # No. 2.: Append the leaf to the left subtree recursively.
            subtrees = self._add_leaf(left, leaf)
            for subtree in subtrees:
                under_threshold, inner_clade = self._create_inner_node(subtree, right, threshold=threshold)
                if under_threshold:
                    if final_iteration:
                        threshold = inner_clade.score
                        min_tree = inner_clade
                    new_trees.append(inner_clade)

            # No. 3.: Append the leaf to the right subtree recursively.
            subtrees = self._add_leaf(right, leaf)
            for subtree in subtrees:
                under_threshold, inner_clade = self._create_inner_node(left, subtree, threshold=threshold)
                if under_threshold:
                    if final_iteration:
                        threshold = inner_clade.score
                        min_tree = inner_clade
                    new_trees.append(inner_clade)

        ####################
        # Extensive search #
        ####################
        else:
            # There are three different kinds of tree formation.

            # No. 1.: create a new root; let the old root be the left subtree and the leaf the right one.
            _, inner_clade = self._create_inner_node(root, leaf, threshold = np.inf)
            new_trees.append(inner_clade)
            if final_iteration and inner_clade.score < threshold:
                threshold = inner_clade.score
                min_tree = inner_clade

            # If the root is a leaf, return current state (it doesn't have subtrees).
            if not root.clades:
                if final_iteration:
                    self.tree = min_tree
                    self.threshold = threshold
                return new_trees

            left = root.clades[0] # Left subtree
            right = root.clades[1] # right subtree

            # No. 2.: Append the leaf to the left subtree recursively.
            subtrees = self._add_leaf(left, leaf)
            for subtree in subtrees:
                _, inner_clade = self._create_inner_node(subtree, right, threshold=np.inf)
                new_trees.append(inner_clade)
                if final_iteration and inner_clade.score < threshold:
                    threshold = inner_clade.score
                    min_tree = inner_clade

            # No. 3.: Append the leaf to the right subtree recursively.
            subtrees = self._add_leaf(right, leaf)
            for subtree in subtrees:
                _, inner_clade = self._create_inner_node(left, subtree, threshold=np.inf)
                new_trees.append(inner_clade)
                if final_iteration and inner_clade.score < threshold:
                    threshold = inner_clade.score
                    min_tree = inner_clade
            

        #####################
        # Change parameters #
        #####################
        # If this is the final iteration and we have found a new best tree, change the appripriate parameters.
        if final_iteration and min_tree is not None:
            self.tree = min_tree
            self.threshold = threshold

        return new_trees

    def _create_inner_node(self, left: ParsimonyClade, right: ParsimonyClade, threshold: int, name: str = None):
        """
        Create new root with given left and right subtree, calculate the parsimony score and create the corresponding sets
        according to Fitch's algorithm. Return prematurely if the score exceeds the given threshold.

        left: root of left subtree
        right: root of right subtree
        threshold: the best result so far (the lower the better)
        name: name of new node (optional)

        returns: whether the resulting parsimony score is under the given threshold; root of new tree (if first parameter is true)
        """

        # Calculate the initial parsimony score by adding the scores of the left and the right node
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
        
    def get_label(self, clade):
        """
        Prettify tree labels.

        clade: current clade
        returns: the label belonging to the current clade
        """
        if clade.name is None:
            return ""
        return clade.name.replace("_", " ")

    def draw_tree(self, show_branch_labels=False):
        """
        Draws tree.

        show_branch_labels: show branch lengths; default is False
        """
        if self.tree is None:
            print("Please first build the tree.")
        else:
            if show_branch_labels:
                draw(self.tree, label_func=self.get_label, branch_labels=(lambda clade: clade.branch_length))
            else:
                draw(self.tree, label_func=self.get_label)
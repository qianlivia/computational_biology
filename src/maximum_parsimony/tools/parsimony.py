import numpy as np
from Bio.Phylo import BaseTree
from Bio.Phylo import draw
from Bio.Align import MultipleSeqAlignment
from collections import Counter
from tools.utils import remove_unprocessed

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

        # Create the dictionary (and list) of taxa
        self.leaves = [] # Leaves as clades
        self.taxa = {} # Node-sequence pairs
        for sequence in alignment:
            leaf = BaseTree.Clade(None, sequence.id)
            self.leaves.append(leaf)
            self.taxa[leaf] = sequence.seq

        self.size = len(self.leaves) # Number of taxa
        self.length = alignment.get_alignment_length() # Length of taxa
        self.bnb = bnb # Branch and bound. If false, use exhaustive search.

    def run(self):
        """
        Run the algorithm.
        """
        # Generate all the trees
        trees = self.generate_trees()

        # Maximum search
        maxi = np.inf
        max_tree = None
        for tree in trees:
            print(trees)
            print()
        self.tree = max_tree

    def generate_trees(self):
        """
        Generate all possible trees.
        """
        trees = [self.leaves[0]]
        if self.size == 1:
            return trees

        for i in range(1, 5):
            curr_leaf = self.leaves[i]
            new_trees = []
            for tree in trees:
                new_trees.extend(self.add_leaf(tree, curr_leaf))
            trees = new_trees

        # TODO optimize
        trees = [self.create_tree(tree) for tree in trees]
        return trees

    def add_leaf(self, root: BaseTree.Clade, leaf: BaseTree.Clade):
        """
        Add leaf to a tree.

        root: root defining the given tree
        leaf: leaf to be added

        returns: trees after the leaf has been added in every possible combination (represented by their roots)
        """

        # There are three different kinds of tree formation.

        # No. 1.: create a new root; let the old root be the left subtree and the leaf the right one.
        # TODO: inner node names?
        inner_clade = self.create_inner_node(root, leaf)
        new_trees = [inner_clade]

        # If the root is a leaf, return current state.
        if not root.clades:
            return new_trees

        left = root.clades[0] # Left subtree
        right = root.clades[1] # right subtree

        # No. 2.: Append the leaf to the left subtree recursively.
        subtrees = self.add_leaf(left, leaf)
        for subtree in subtrees:
            inner_clade = self.create_inner_node(subtree, right)
            new_trees.append(inner_clade)

        # No. 3.: Append the leaf to the right subtree recursively.
        subtrees = self.add_leaf(right, leaf)
        for subtree in subtrees:
            inner_clade = self.create_inner_node(left, subtree)
            new_trees.append(inner_clade)

        return new_trees

    def create_tree(self, root: BaseTree.Clade):
        """
        Create tree with given root.

        root: root

        returns: tree
        """
        return BaseTree.Tree(root, rooted=False)

    def create_inner_node(self, left: BaseTree.Clade, right: BaseTree.Clade, name: str = None):
        """
        Create root with given left and right subtree.

        left: root of left subtree
        right: root of right subtree

        returns: root of new tree
        """
        root = BaseTree.Clade(None, name)
        root.clades.append(left)
        root.clades.append(right)
        return root

    def calc_parsimony(self, tree: BaseTree.Tree, parents: dict):
        """
        Fitch algorithm. Calculate score for every site (self.length).
        
        tree: tree
        parents: dictionary of parent nodes
        """

        sum = 0
        for u in range(self.length):
            sum += self.parsimony_for_site(tree, parents, u)

        return sum

    def parsimony_for_site(self, tree: BaseTree, parents: dict, u: int):
        """
        Calculate parsimony score for given site. From bottom to top.
        https://www.cs.helsinki.fi/bioinformatiikka/mbi/courses/07-08/itb/slides/itb0708_slides_158-191.pdf
        https://www.bio.fsu.edu/~stevet/BSC5936/Swofford.F2003.2.pdf
        https://tel.archives-ouvertes.fr/tel-01479049/document
        http://pages.stat.wisc.edu/~larget/Genetics629/Fall2009/outline2.pdf

        tree: tree
        parents: dictionary of parent nodes
        u: given site
        """
        score = 0 # Parsimony score.
        sets_of_bases = {} # Dictionary of sets. Format: {Leaf0: {'A'}, Leaf1: {'T'}, Leaf2: {'C'}, ..., InnerNode0: {'A', 'T', 'G'}}

        # Create the initial sets for leaves.
        nodes_to_be_processed = []
        for leaf in self.leaves:
            base = self.taxa[leaf][u]
            sets_of_bases[leaf] = {base}
            parent = parents[leaf]
            nodes_to_be_processed.append(parent)

        # Remove clades that do not occur as many times as how many children they have.
        nodes_to_be_processed = remove_unprocessed(nodes_to_be_processed)

        # Do this until we reach the root.
        while nodes_to_be_processed:
            next_nodes = []
            # Iterate over inner nodes.
            for node in nodes_to_be_processed:
                children = node.clades
                X = set() # Intersection
                Y = set() # Union

                # Find the intersection and the union of the children sets.
                for child in children:
                    set_of_base = sets_of_bases[child]
                    X = X.intersection(set_of_base)
                    Y = Y.union(set_of_base)
                # If not empty, X is the set belonging to the current node.
                if X:
                    sets_of_bases[node] = X
                # Otherwise, the set for the node is the union of the children sets and one is added to the score.
                else:
                    sets_of_bases[node] = Y
                    score += 1

                # Check if current node has a parent.
                if node in parents:
                    # Add the parent of the current node to the list of nodes to be processed.
                    next_nodes.append(parent)

            # Remove clades that do not occur as many times as how many children they have.
            nodes_to_be_processed = remove_unprocessed(next_nodes)

        return score
        
    def get_label(self, clade):
        """
        Prettify tree labels.

        clade: current clade
        returns: the label belonging to the current clade
        """
        if clade.name.startswith("Inner"):
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
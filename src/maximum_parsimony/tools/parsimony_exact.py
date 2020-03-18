import numpy as np
from Bio.Align import MultipleSeqAlignment
from tools.parsimony import Parsimony, ParsimonyClade

class ParsimonyExact(Parsimony):
    """
    Class for exact methods using maximum parsimony (branch and bound, exhaustive search).
    """

    def __init__(self, alignment: MultipleSeqAlignment, bnb: bool = True):
        """
        Constructor.

        alignment: MultipleSeqAlignment containing the alignment
        bnb: decides whether to use branch andbound (optional). Default is True. 
        """
        super(ParsimonyExact, self).__init__(alignment)
        self.bnb = bnb # Branch and bound. If false, use exhaustive search.

    def run(self, print_best: bool = False):
        """
        Generate all possible trees and find the one with the lowest score.

        print_best: print best tree (optional). Default is False.
        """

        super(ParsimonyExact, self).run()

        if self.size <= 0:
            print("There aren't enough taxa.")
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
        """
        for tree in self.trees:
            print(tree)
            print()
        """

        print("Number of trees considered as local maxima:", len(self.trees))

        if print_best:
            print(self.tree)

    def _DFS_trees(self, queue: list, leaf_no: int):
        """
        Investigate every possible tree using depth-first search. Add new taxa to the current trees one by one in every possible combination and return.
        Store best trees if in the last iteration.

        queue: queue containing the trees to be extended
        leaf_no: the index of the leaf to extend the trees with
        """
        while queue:
            tree = queue.pop(0)

            if leaf_no < self.size - 1:
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
        threshold: the best result so far; the lower the better (optional). Default is infinity.
        final_iteration: decides whether this is the final iteration when the last taxon is added (optional). Default is False.

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
            _, inner_clade = self._create_inner_node(root, leaf, threshold=np.inf)
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
        # Update parameters #
        #####################
        # If this is the final iteration and we have found a new best tree, change the appripriate parameters.
        if final_iteration and min_tree is not None:
            self.tree = min_tree
            self.threshold = threshold

        return new_trees
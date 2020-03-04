import numpy as np
from Bio.Align import MultipleSeqAlignment
from tools.parsimony import Parsimony, ParsimonyClade

class ParsimonyHeuristics(Parsimony):
    """
    Class for heuristic methods using maximum parsimony (subtree pruning and regrafting).
    """

    def __init__(self, alignment: MultipleSeqAlignment, seed: int = 0, limit: int = np.inf):
        """
        Constructor.

        alignment: MultipleSeqAlignment containing the alignment
        """
        super(ParsimonyHeuristics, self).__init__(alignment)
        np.random.seed(seed)
        self.limit = limit # Maximum number of trees to investigate
        self.parents = {} # Dictionary for storing the parents of nodes
        self.inner_nodes = set() # Set of inner nodes
        self.inner_nodes_wo_root = set() # Set of inner nodes without the root

    def run(self, print_best = False):
        """
        Generate all possible trees and find the one with the lowest score.

        print_best: print best tree (optional). Default is False.
        """

        super(ParsimonyHeuristics, self).run()
        
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
        self._generate_initial_tree()
        # self._SPR()
        self.tree = self.create_tree(self.tree)
        
        # Print tree
        if print_best:
            print(self.tree)

        print(self.inner_nodes)
        print()
        print(self.inner_nodes_wo_root)
        print()
        print(self.parents)

    def _generate_initial_tree(self):
        leaves = self.leaves
        np.random.shuffle(leaves)
        root = leaves[0]
        for i in range(1, self.size):
            # Create new inner node
            _, inner_node = self._create_inner_node(root, leaves[i], name="Inner{}".format(i - 1), threshold=np.inf)

            # Store parents
            self.parents[root] = inner_node
            self.parents[leaves[i]] = inner_node

            # Store inner node
            self.inner_nodes.add(inner_node)
            self.inner_nodes_wo_root.add(inner_node)

            root = inner_node

        # Copy set of inner nodes
        self.inner_nodes_wo_root.discard(root)
        # Remove root from the set of inner nodes (last element)
        
        self.tree = root

    def _SPR(self):
        # TODO: implement exit condition dependent on change in parsimony score
        if self.limit == np.inf:
            self.limit = 1000

        iter = 0
        while iter < self.limit:
            # Choose a random edge that's not connected to the root; first pick a random inner node and then pick one of its children.
            parent = np.random.choice(self.inner_nodes_wo_root)
            child_ind = np.random.choice([0, 1])
            child = parent.clades[child_ind]

            # Take this edge out and rewire accordingly.
            new_parent = self.parents[parent] # grandparent
            new_child = parent.clades[1 - child_ind] # grandchild

            # Create new edges.
            if new_parent.clades[0] == parent:
                new_parent.clades[0] = new_child
            else:
                new_parent.clades[1] = new_child
            self.parents[new_child] = new_parent

            # Delete old edges.
            self.parents[parent] = None
            parent.clades.pop(1 - child_ind)

            # Choose another node and one of its children (excluding subtree we took out and the edge that we just created).
            # First, collect all the inner nodes we shouldn't consider choosing from.
            used_nodes = {parent, child}
            curr_nodes = [child]
            while curr_nodes:
                curr = curr_nodes.pop(0)
                if curr.clades:
                    used_nodes.add(curr.clades[0])
                    used_nodes.add(curr.clades[1])
                    curr_nodes.extend(curr.clades)

            # Then choose from the remainder.
            nodes_to_choose_from = self.inner_nodes - used_nodes
            chosen_parent = np.random.choice(nodes_to_choose_from)

            #If this is the parent of the original edge, choose the other edge belonging to this node:
            if chosen_parent == new_parent:
                pass

            # Insert subtree.
            # Calculate new score.

            iter += 1

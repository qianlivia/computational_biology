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
        np.random.seed(seed) # Set seed
        self.limit = limit # Maximum number of trees to investigate
        self.parents = {} # Dictionary for storing the parents of nodes
        self.inner_nodes = set() # Set of inner nodes
        self.inner_nodes_wo_root = set() # Set of inner nodes without the root

    def run(self, print_best: bool = False):
        """
        Generate all possible trees and find the one with the lowest score.

        print_best: print best tree (optional). Default is False.
        """

        super(ParsimonyHeuristics, self).run()
        
        if self.size <= 0:
            print("There aren't enough taxon.")
            return

        # Return if there is only one taxon
        if self.size == 1:
            inner_clade = ParsimonyClade(None, None, sets=[], score=0)
            inner_clade.clades.append(self.leaves[0])
            self.tree = self.create_tree(inner_clade)
            return
            
        # Initialize
        self._generate_initial_tree()

        # Return if there are only two taxa
        if self.size == 2:
            self.tree = self.create_tree(self.tree)
            return

        # Run SPR
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

        # Remove root from the set of inner nodes (last element)
        self.inner_nodes_wo_root.discard(root)
        
        self.tree = root

    def _SPR(self):
        """
        Subtree pruning and regrafting. Prune and reallocate a random subtree for a certain number of iterations (self.limit).
        Variables:

        parent: parent of the root of the subtree that needs to be moved - will be moved together with the subtree
        child: root of the subtree that needs to be moved
        new_parent: grandparent of the root of the subtree that needs to be moved
        new_child: the other child of variable "parent" - this will replace "parent" in the hierarchy
        chosen_parent: will be the new parent of variable "parent" 
        chosen_child: will be the new child of variable "parent", along with the subtree that we pruned and regrafted

        """
        # TODO: implement exit condition dependent on change in parsimony score
        if self.limit == np.inf:
            self.limit = 100

        iter = 0
        while iter < self.limit:
            # Choose a random edge that's not connected to the root; first pick a random inner node, then pick one of its children.
            parent = np.random.choice(self.inner_nodes_wo_root)
            child_ind = np.random.choice([0, 1])
            child = parent.clades[child_ind]

            # Take this edge out and rewire accordingly.
            new_parent = self.parents[parent] # grandparent
            new_child = parent.clades[1 - child_ind] # grandchild

            # Create new edges.
            is_left = True
            if new_parent.clades[0] == parent:
                new_parent.clades[0] = new_child
            else:
                is_left = False
                new_parent.clades[1] = new_child
            self.parents[new_child] = new_parent

            # Delete old edges.
            self.parents[parent] = None
            parent.clades.pop(1 - child_ind)

            # Choose another node and one of its children (excluding the subtree we took out and the edge that we just created).
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

            # If this is the parent of the original edge, choose the other edge belonging to this node:
            if chosen_parent == new_parent:
                if is_left:
                    chosen_child_ind = 1
                else:
                    chosen_child_ind = 0
            else:
                chosen_child = np.random.choice([0, 1])
                chosen_child = chosen_parent.clades[chosen_child_ind]

            # Insert subtree.
            chosen_parent.clades[chosen_child_ind] = parent
            self.parents[parent] = chosen_parent
            if chosen_child_ind: # if the chosen child is on the right
                parent.clades.append(chosen_child)
            else: # if it is on the left
                parent.clades.insert(0, chosen_child)
            self.parents[chosen_child] = parent

            # Calculate new parsimony score.
            self._recalculate_score(new_parent, parent, self.tree)

            iter += 1

    def _recalculate_score(self, first_node: ParsimonyClade, second_node: ParsimonyClade, root: ParsimonyClade):
        """
        Calculate new score by detecting the nodes that need to be changed because of the regrafting.
        Trace the first node back to the root and put its ancestors into a list.
        Trace the second one; check whether its ancestors are already in the list created previously. The first ancestor
        to be present in that list is the first intersection of the two paths; this means that this and all the successive
        nodes need to be processed after the nodes that are only present in one of the two paths.

        first_node: the first node of the first path
        second_mode: the first node of the second path
        root: the root of the tree
        """

        # Create first path
        first_path = [first_node]
        curr = new_parent
        while curr != root:
            curr = self.parents[curr]
            first_path.append(curr)

        # Create second path
        second_path = []
        curr = second_node
        
        # If first element exists in first_path
        if curr in first_path:
            ind = first_path.index(curr)
            common_path = first_path[ind:]
            first_path = first_path[:ind]
        # If first element is not in first_path
        else:
            # Initialize second_path with the first element
            second_path.append(curr)

            while curr != root:
                curr = self.parents[curr]
                if curr in first_path:
                    ind = first_path.index(curr)
                    common_path = first_path[ind:]
                    first_path = first_path[:ind]
                    break
                else:
                    second_path.append(curr)

        # Change the parsimony score in the elements of all three lists
        self._calculate_score_for_list(first_path)
        self._calculate_score_for_list(second_path)
        self._calculate_score_for_list(common_path)
        
    def _calculate_score_for_list(self, lst: list):
        """
        Iterate over a list and change the parsimony score and sets of all the elements in the list.

        lst: list
        """

        for elem in lst:
            # Initialize variables
            sets = []
            score = 0
            left = elem.clades[0]
            right = elem.clades[1]
            
            # Count parsimony score for every site (self.length elements)
            for u in range(self.length):
                u_score, u_set = self._parsimony_for_site(left, right, u)
                score += u_score
                sets.append(u_set)

            elem.score = score
            elem.sets = sets
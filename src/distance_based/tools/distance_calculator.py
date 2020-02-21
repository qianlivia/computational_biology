from typing import List
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo import BaseTree
from Bio.Phylo import draw
import itertools
from tools.utils import Matrix, Vector

class Distance_Calculator():
    """
    Class containing distance-based methods.
    """

    def __init__(self, mode='NJ'):
        """ Constructor.

        mode: NJ, UPGMA or WPGMA
        """
        self.mode = mode
        self.tree = None

    def _pairwise_distance(self, sequence1, sequence2):
        """
        Determines the pairwise distance between two sequences.

        sequence1: first sequence
        sequence2: second sequence
        returns: the pairwise distance
        """
        sum = 0
        pairs = zip(sequence1, sequence2)
        length = 0
        for pair in pairs:
            if (pair[0] != pair[1]):
                sum += 1
            length += 1
        return sum

    def create_distance_matrix(self, alignment: MultipleSeqAlignment):
        """
        Create distance matrix by looking at pairwise differences.

        alignment: MultipleSeqAlignment containing the alignment
        returns: distance matrix
        """
        names = [sequence.id for sequence in alignment]
        dm = Matrix(names)
        for sequence1, sequence2 in itertools.combinations(alignment, 2):
            dm[sequence1.id][sequence2.id] = self._pairwise_distance(sequence1, sequence2)
        return dm

    def build_tree(self, distance_matrix: Matrix):
        """
        Build tree.

        distance_matrix: distance matrix
        """
        dm = distance_matrix.copy()
        if self.mode == 'NJ':
            self._neigbor_joining(dm)
        if self.mode == 'UPGMA':
            self._upgma(dm)
        if self.mode == 'WPGMA':
            self._wpgma(dm)

    def _neigbor_joining(self, dm: Matrix):
        """ Neighbor joining.

            dm: distance matrix
        """
        # Create tree nodes
        clades = [BaseTree.Clade(None, name) for name in dm.labels]
        nr_inner = 0

        N = dm.size

        # Run until two nodes remain.
        while N > 2:
            labels = dm.labels
            # Divergence
            divergence = Vector(labels)
            for label in labels:
                # Sum up all the distances for the given label
                divergence[label] = sum([dm[label][i] for i in labels]) + sum([dm[i][label] for i in labels])

            # New distance matrix
            dm_new = Matrix(labels)
            for label1, label2 in itertools.combinations(labels, 2):
                dm_new[label1][label2] = dm[label1][label2] - (divergence[label1] + divergence[label2]) / (N - 2)

            # Find the index of the smallest value in the distance matrix
            argmin1, argmin2 = dm_new.argmin() # Labels belonging to the smallest value
            minpos1, minpos2 = dm_new.posmin() # Indices belonging to the smallest value
            c1, c2 = clades[minpos1], clades[minpos2] # Fetch the corresponding clades

            # New inner node
            node_name = "Inner{}".format(nr_inner)
            inner_clade = BaseTree.Clade(None, node_name)

            # Calculate branch length for the two old nodes
            c1.branch_length = dm[argmin1][argmin2]/2 + (divergence[argmin1] - divergence[argmin2]) / (2 * (N - 2))
            c2.branch_length = dm[argmin1][argmin2] - c1.branch_length

            # Append these to the new node
            inner_clade.clades.append(c1)
            inner_clade.clades.append(c2)

            # Calculate the distance from the new node to all the other nodes (not including the ones we want to remove).
            dm.add([node_name])
            
            neighbor_labels = set(labels) - {node_name, argmin1, argmin2}
            for label in neighbor_labels:
                dm[label][node_name] = (dm[argmin1][label] + dm[label][argmin1] + dm[argmin2][label] + dm[label][argmin2] - dm[argmin1][argmin2]) / 2

            # Delete appropriate rows, columns and labels.
            dm.drop([argmin1, argmin2])

            # Discard the old nodes and append the new one
            del clades[minpos2]
            del clades[minpos1]
            clades.append(inner_clade)
            nr_inner += 1

            # Number of nodes remaining
            N = dm.size

        # Create last inner node
        node_name = "Inner{}".format(nr_inner)
        inner_clade = BaseTree.Clade(None, node_name)

        # Append last two nodes
        inner_clade.clades.append(clades[0])
        inner_clade.clades.append(clades[1])

        # Create tree
        self.tree = BaseTree.Tree(inner_clade, rooted=False)

    def _wpgma(self, dm: Matrix):
        """ WPGMA.

            dm: distance matrix
        """
        # Create tree nodes
        clades = [BaseTree.Clade(None, name) for name in dm.labels]
        clade_children_branch_length = Vector(dm.labels)
        nr_inner = 0

        # Number of nodes
        N = dm.size
        
        # Run until one node remains.
        while N > 1:
            # Find the index of the smallest value in the distance matrix
            argmin1, argmin2 = dm.argmin(only_positive=True) # Labels belonging to the smallest value
            minpos1, minpos2 = dm.posmin(only_positive=True) # Indices belonging to the smallest value
            c1, c2 = clades[minpos1], clades[minpos2] # Fetch the corresponding clades

            # New inner node
            node_name = "Inner{}".format(nr_inner)
            inner_clade = BaseTree.Clade(0, node_name)

            # Calculate branch length for the two old nodes
            c1.branch_length = dm[argmin1][argmin2]/2 - clade_children_branch_length[argmin1]
            c2.branch_length = dm[argmin1][argmin2]/2 - clade_children_branch_length[argmin2]

            # Append these to the new node
            inner_clade.clades.append(c1)
            inner_clade.clades.append(c2)

            # Append new node to the appropriate places
            dm.add([node_name])
            clade_children_branch_length.add([node_name])
            clade_children_branch_length[node_name] = dm[argmin1][argmin2]/2

            # Calculate the distance from the new node to all the other nodes (not including the ones we want to remove).
            neighbor_labels = set(dm.labels) - {node_name, argmin1, argmin2}
            for label in neighbor_labels:
                # Either dm[argmin1][label] or dm[label][argmin1] is zero. Same with dm[argmin2][label] and dm[label][argmin2].
                dm[label][node_name] = (dm[argmin1][label] + dm[label][argmin1] + dm[argmin2][label] + dm[label][argmin2]) / 2

            # Delete appropriate rows, columns and labels.
            dm.drop([argmin1, argmin2])

            # Discard the old nodes and append the new one
            del clades[minpos2]
            del clades[minpos1]
            clades.append(inner_clade)
            nr_inner += 1

            # Number of nodes remaining
            N = dm.size

        # Create tree
        self.tree = BaseTree.Tree(clades[0], rooted=False)

    def _upgma(self, dm: Matrix):
        """ UPGMA.

            dm: distance matrix
        """
        # Create tree nodes
        clades = [BaseTree.Clade(None, name) for name in dm.labels]
        clade_children_branch_length = Vector(dm.labels)
        clade_nr_of_children = Vector(dm.labels, init_value=1)
        nr_inner = 0

        # Number of nodes
        N = dm.size
        
        # Run until one node remains.
        while N > 1:
            # Find the index of the smallest value in the distance matrix
            argmin1, argmin2 = dm.argmin(only_positive=True) # Labels belonging to the smallest value
            minpos1, minpos2 = dm.posmin(only_positive=True) # Indices belonging to the smallest value
            c1, c2 = clades[minpos1], clades[minpos2] # Fetch the corresponding clades

            # New inner node
            node_name = "Inner{}".format(nr_inner)
            inner_clade = BaseTree.Clade(0, node_name)

            # Calculate branch length for the two old nodes
            c1.branch_length = dm[argmin1][argmin2]/2 - clade_children_branch_length[argmin1]
            c2.branch_length = dm[argmin1][argmin2]/2 - clade_children_branch_length[argmin2]

            # Append these to the new node
            inner_clade.clades.append(c1)
            inner_clade.clades.append(c2)

            # Append new node
            dm.add([node_name])
            
            # Branch length from the current node to the leaves
            clade_children_branch_length.add([node_name])
            clade_children_branch_length[node_name] = dm[argmin1][argmin2]/2

            # Number of children in the new group
            clade_nr_of_children.add([node_name])
            clade_nr_of_children[node_name] = clade_nr_of_children[argmin1] + clade_nr_of_children[argmin2]

            # Calculate the distance from the new node to all the other nodes (not including the ones we want to remove).
            neighbor_labels = set(dm.labels) - {node_name, argmin1, argmin2}
            for label in neighbor_labels:
                # Either dm[argmin1][label] or dm[label][argmin1] is zero. Same with dm[argmin2][label] and dm[label][argmin2].
                w1, w2 = clade_nr_of_children[argmin1], clade_nr_of_children[argmin2]
                dm[label][node_name] = \
                    ((dm[argmin1][label] + dm[label][argmin1]) * w1 + (dm[argmin2][label] + dm[label][argmin2]) * w2) / (w1 + w2)

            # Delete appropriate rows, columns and labels.
            dm.drop([argmin1, argmin2])
            clade_children_branch_length.drop([])

            # Discard the old nodes and append the new one. The old clades need to be discarded because relative position is used when referring to them!
            del clades[minpos2]
            del clades[minpos1]
            clades.append(inner_clade)
            nr_inner += 1

            # Number of nodes remaining
            N = dm.size

        # Create tree
        self.tree = BaseTree.Tree(clades[0], rooted=False)

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
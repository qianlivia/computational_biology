from typing import List
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped
from Bio.Phylo import BaseTree
from Bio.Phylo import draw
import itertools
from utils import Matrix, Vector

class Distance_Calculator():

    def __init__(self, mode='NJ'):
        """ Mode: NJ, UPGMA, WPGMA. """
        self.mode = mode
        self.tree = None

    def _pairwise_distance(self, sequence1, sequence2):
        sum = 0
        pairs = zip(sequence1, sequence2)
        length = 0
        for pair in pairs:
            if (pair[0] != pair[1]):
                sum += 1
            length += 1
        return sum

    def create_distance_matrix(self, alignment: MultipleSeqAlignment):
        names = [sequence.id for sequence in alignment]
        dm = Matrix(names)
        for sequence1, sequence2 in itertools.combinations(alignment, 2):
            dm[sequence1.id][sequence2.id] = self._pairwise_distance(sequence1, sequence2)
        return dm

    def build_tree(self, distance_matrix: Matrix):
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
            dm.add(node_name)
            
            neighbor_labels = set(labels) - {node_name, argmin1, argmin2}
            for label in neighbor_labels:
                dm[label][node_name] = (dm[argmin1][label] + dm[label][argmin1] + dm[argmin2][label] + dm[label][argmin2] - dm[argmin1][argmin2]) / 2

            # Delete appropriate rows, columns and labels.
            dm.drop([argmin1, argmin2])

            # Replace one of the nodes with the new node and discard the other one
            clades[minpos1] = inner_clade
            del clades[minpos2]
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

    def get_label(self, leaf):
        """
        Prettify tree labels.
        leaf: current leaf
        """
        if leaf.name.startswith("Inner"):
            return ""
        return leaf.name.replace("_", " ")

    def _upgma(self, dm: Matrix):
        """ UPGMA.
            dm: distance matrix
        """
        # Create tree nodes
        clades = [BaseTree.Clade(0, name) for name in dm.labels]
        nr_inner = 0

        # Number of nodes
        N = dm.size
        
        # Run until two nodes remain.
        while N > 2:
            # Find the index of the smallest value in the distance matrix
            argmin1, argmin2 = dm.argmin(exclude_zeros=True) # Labels belonging to the smallest value
            minpos1, minpos2 = dm.posmin(exclude_zeros=True) # Indices belonging to the smallest value
            c1, c2 = clades[minpos1], clades[minpos2] # Fetch the corresponding clades

            # New inner node
            node_name = "Inner{}".format(nr_inner)
            inner_clade = BaseTree.Clade(0, node_name)

            # Calculate branch length for the two old nodes
            # TODO this needs to be rewritten
            c1.branch_length = dm[argmin1][argmin2]/2 - clades[minpos1].branch_length
            c2.branch_length = dm[argmin1][argmin2]/2 - clades[minpos2].branch_length
            

            # Append these to the new node
            inner_clade.clades.append(c1)
            inner_clade.clades.append(c2)

            # Calculate the distance from the new node to all the other nodes (not including the ones we want to remove).
            dm.add(node_name)

            neighbor_labels = set(dm.labels) - {node_name, argmin1, argmin2}
            for label in neighbor_labels:
            # TODO this needs to be rewritten
                dm[label][node_name] = (dm[argmin1][label] + dm[label][argmin1] + dm[argmin2][label] + dm[label][argmin2]) / 2

            # Delete appropriate rows, columns and labels.
            dm.drop([argmin1, argmin2])

            # Replace one of the nodes with the new node and discard the other one
            clades[minpos1] = inner_clade
            del clades[minpos2]
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
        pass

    def draw_tree(self):
        if self.tree is None:
            print("Please first build the tree.")
        else:
            draw(self.tree, label_func=self.get_label)
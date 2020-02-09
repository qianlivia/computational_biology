from typing import List
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped
from Bio.Phylo import BaseTree
from tabulate import tabulate
import itertools
import numpy as np
import pandas as pd

class Vector():
    
    def __init__(self, labels: str):
        self.labels = labels
        self.size = len(labels)
        self.data = pd.DataFrame(np.zeros((self.size)))
        self.data.index = labels

    def __getitem__(self, item):
        return self.data[0][item]

    def __setitem__(self, item, value):
        self.data[0][item] = value

    def __str__(self):
        return tabulate(self.data, headers='keys', tablefmt='psql')

    def argmax(self):
        return self.data.max(axis=1).idxmax()

    def argmin(self):
        return self.data.min(axis=1).idxmin()

class Matrix():

    def __init__(self, labels: str):
        self.labels = labels
        self.size = len(labels)
        self.data = pd.DataFrame(data=np.zeros((self.size, self.size)), index=labels, columns=labels)

    def __getitem__(self, item):
        return self.data.__getitem__(item)

    def __setitem__(self, item, value):
        self.data.__setitem__(item, value)

    def __str__(self):
        return tabulate(self.data, headers='keys', tablefmt='psql')

    def argmax(self):
        return (self.data.max(axis=0).idxmax(), self.data.max(axis=1).idxmax())

    def argmin(self):
        return (self.data.min(axis=0).idxmin(), self.data.min(axis=1).idxmin())

    def posmax(self):
        a, b = np.unravel_index(np.argmax(self.data.values, axis=None), self.data.values.shape)
        return b, a

    def posmin(self):
        a, b = np.unravel_index(np.argmin(self.data.values, axis=None), self.data.values.shape)
        return b, a

class Distance_Calculator():

    def __init__(self, mode='NJ'):
        """ Mode: NJ, UPGMA, WPGMA. """
        self.mode = mode

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
        # Create tree nodes
        clades = [BaseTree.Clade(None, name) for name in dm.labels]
        nr_inner = 0

        # Divergence
        N = dm.size

        while N > 0:
            divergence = Vector(dm.labels)
            for label in dm.labels:
                divergence[label] = sum([dm[label][i] for i in dm.labels]) + sum([dm[i][label] for i in dm.labels])
            print(divergence)

            # New distance matrix
            dm_new = Matrix(dm.labels)
            for label1, label2 in itertools.combinations(dm.labels, 2):
                dm_new[label1][label2] = dm[label1][label2] - (divergence[label1] + divergence[label2]) / (N - 2)

            # Find the index of the smallest value in the distance matrix
            argmin1, argmin2 = dm_new.argmin()
            minpos1, minpos2 = dm_new.posmin()
            c1, c2 = clade[minpos1], clade[minpos2]

            # New inner node
            node_name = "Inner{}".format(nr_inner)
            inner_clade = BaseTree.Clade(None, node_name)

            # Calculate branch length for the two old nodes
            c1.branch_length = dm[argmin1][argmin2]/2 + (divergence[argmin1] + divergence[argmin2]) / (2 * (N - 2))
            c2.branch_length = dm[argmin1][argmin2] - c1.branch_length

            # Append these to the new node
            inner_clade.clades.append(c1)
            inner_clade.clades.append(c2)

            # Calculate the distance from the new node to all the other nodes
            dm.labels.remove(argmin1)
            dm.labels.remove(argmin2)

            for label in dm.labels:
                dm[node_name][label] = (dm[argmin1][label] + dm[argmin2][label] - dm[argmin1][argmin2]) / 2

            del dm[argmin1]
            del dm[argmin2]
            # TODO delete rows, too!

            # Replace one of the nodes with the new node and discard the other one
            clades[minpos1] = inner_clade
            del clades[minpos2]
            nr_inner += 1

            dm.labels.append(node_name)
            N = dm.size

    def _upgma(self, dm: Matrix):
        pass

    def _wpgma(self, dm: Matrix):
        pass
from typing import List
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped
from tabulate import tabulate
import itertools
import numpy as np
import pandas as pd

class Vector():
    
    def __init__(self, labels: str):
        self.labels = labels
        self.size = len(labels)
        self._vector = np.zeros((self.size))
        self.data = pd.DataFrame(self._vector)
        self.data.index = labels

    def __getitem__(self, item):
        return self.data[0][item]

    def __setitem__(self, item, value):
        self.data[0][item] = value

    def __str__(self):
        return tabulate(self.data, headers='keys', tablefmt='psql')

class Matrix():

    def __init__(self, labels: str):
        self.labels = labels
        self.size = len(labels)
        self._matrix = np.zeros((self.size, self.size))
        self.data = pd.DataFrame(data=self._matrix, index=labels, columns=labels)

    def __getitem__(self, item):
        return self.data.__getitem__(item)

    def __setitem__(self, item, value):
        self.data.__setitem__(item, value)

    def __str__(self):
        return tabulate(self.data, headers='keys', tablefmt='psql')

class Distance_Calculator():

    def __init__(self, mode='NJ'):
        """ Mode: NJ, UPGMA, WPGMA. """
        self.mode = mode

    def _pairwise_distance(self, sequence1, sequence2):
        sum = 0
        pairs = zip(sequence1, sequence2)
        length = 0
        for pair in pairs:
            if (pair[0] == pair[1]):
                sum += 1
            length += 1
        return 1 - (sum/length)

    def create_distance_matrix(self, alignment: MultipleSeqAlignment):
        names = [sequence.id for sequence in alignment]
        dm = Matrix(names)
        for sequence1, sequence2 in itertools.combinations(alignment, 2):
            dm[sequence1.id, sequence2.id] = self._pairwise_distance(sequence1, sequence2)
        return dm

    def build_tree(self, dm: Matrix):
        if self.mode == 'NJ':
            # Divergence
            divergence = Vector(dm.labels)
            for label in dm.labels:
                print(label, sum([dm[label][i] for i in dm.labels]) + sum([dm[i][label] for i in dm.labels]))
                divergence[label] = sum([dm[label][i] for i in dm.labels]) + sum([dm[i][label] for i in dm.labels])
            print(divergence)
            # New distance matrix
            dm_new = Matrix(dm.labels)
            for label1, label2 in itertools.combinations(dm.labels, 2):
                dm_new[label1][label2] = dm[label1][label2] - [divergence[label1] + divergence[label2]] / (dm.size - 2)
            print(dm_new)
        if self.mode == 'UPGMA':
            pass

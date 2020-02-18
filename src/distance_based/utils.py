from typing import List
import numpy as np
import pandas as pd
from tabulate import tabulate

class Vector():
    
    def __init__(self, labels: str):
        pd.options.mode.chained_assignment = None
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
        pd.options.mode.chained_assignment = None
        self.labels = labels
        self.size = len(labels)
        self.data = pd.DataFrame(data=np.zeros((self.size, self.size)), index=labels, columns=labels)

    def __getitem__(self, item):
        return self.data.__getitem__(item)

    def __setitem__(self, item, value):
        self.data.__setitem__(item, value)

    def __str__(self):
        return tabulate(self.data, headers='keys', tablefmt='psql')

    def copy(self):
        obj = type(self).__new__(self.__class__)
        obj.__dict__.update(self.__dict__)
        return obj

    def add(self, label):
        # Insert column.
        self.data[label] = 0

        # Add label to the list of labels and increase size.
        self.labels.append(label)
        self.size = len(self.labels)

        # Insert row.
        self.data.loc[label] = [0] * (self.size)

    def drop(self, labels):
        self.data.drop(labels, inplace=True, axis=0)
        self.data.drop(labels, inplace=True, axis=1)
        for l in labels:
            self.labels.remove(l)
        self.size = len(self.labels)

    def argmax(self):
        return (self.data.max(axis=0).idxmax(), self.data.max(axis=1).idxmax())

    def argmin(self, exclude_zeros=False):
        if exclude_zeros:
            mask = (self.data == 0)
            masked_data = self.data.mask(mask, np.inf, inplace=False)
            return (masked_data.min(axis=0).idxmin(), masked_data.min(axis=1).idxmin())
        return (self.data.min(axis=0).idxmin(), self.data.min(axis=1).idxmin())

    def posmax(self):
        a, b = np.unravel_index(np.argmax(self.data.values, axis=None), self.data.values.shape)
        return b, a

    def posmin(self, exclude_zeros=False):
        if exclude_zeros:
            mask = (self.data == 0)
            masked_data = self.data.mask(mask, np.inf, inplace=False)
            a, b = np.unravel_index(np.argmin(masked_data.values, axis=None), masked_data.values.shape)
        else:
            a, b = np.unravel_index(np.argmin(self.data.values, axis=None), self.data.values.shape)
        return b, a
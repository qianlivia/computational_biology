from typing import List
import numpy as np
import pandas as pd
from tabulate import tabulate

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

    def drop_rows(self, labels):
        self.data.drop(labels, inplace=True, axis=0)
        
    def drop_columns(self, labels):
        self.data.drop(labels, inplace=True, axis=1)

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
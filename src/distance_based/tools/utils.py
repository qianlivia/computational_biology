from typing import List
import numpy as np
import pandas as pd
from tabulate import tabulate

class Vector():
    """
    Vector capable of handling strings as indices.
    """
    
    def __init__(self, labels: List, init_value: float = 0.0):
        """
        Constructor.

        labels: labels for the resulting frame
        init_value: default value of elements
        """
        pd.options.mode.chained_assignment = None
        self.labels = labels.copy()
        self.size = len(self.labels)
        self.data = pd.DataFrame(init_value * np.ones((self.size)))
        self.data.index = self.labels
        self.init_value = init_value

    def __getitem__(self, item: int):
        """
        Returns the element at a certain position.

        item: position

        returns: the value of the chosen item
        """
        return self.data[0][item]

    def __setitem__(self, item: int, value: float):
        """
        Sets the element at a certain position to a certain value.

        item: position
        value: new value
        """
        self.data[0][item] = value

    def __str__(self):
        """
        The string form of the vector. Contains labels.

        returns: the prettified string
        """
        return tabulate(self.data, headers='keys', tablefmt='psql')
        
    def copy(self):
        """
        Copy current instance.

        returns: new instance
        """
        obj = type(self).__new__(self.__class__)
        obj.__dict__.update(self.__dict__)
        return obj

    def add(self, labels: List):
        """
        Add new labels and corresponding rows.

        labels: labels to add
        """
        for label in labels:
           # Add labels to the list of labels and increase size.
            self.labels.append(label)

            # Insert row.
            self.data.loc[label] = [self.init_value]
            
        self.size = len(self.labels)

    def drop(self, labels: List):
        """
        Drop labels and corresponding rows.

        labels: labels to drop
        """
        self.data.drop(labels, inplace=True, axis=0)
        for l in labels:
            self.labels.remove(l)
        self.size = len(self.labels)

    def argmax(self):
        """
        Returns the label belonging to the maximum value.
        """
        return self.data.max(axis=1).idxmax()

    def argmin(self):
        """
        Returns the label belonging to the minimum value.
        """
        return self.data.min(axis=1).idxmin()

class Matrix():
    """
    Matrix capable of handling strings as indices.
    """

    def __init__(self, labels: str, init_value: float = 0.0):
        """
        Constructor.

        labels: labels for the resulting frame
        init_value: default value of elements
        """
        pd.options.mode.chained_assignment = None
        self.labels = labels.copy()
        self.size = len(self.labels)
        self.data = pd.DataFrame(data=(init_value * np.ones((self.size, self.size))), index=labels, columns=labels)
        self.init_value = init_value

    def __getitem__(self, item: int):
        """
        Returns the element at a certain position.

        item: position

        returns: the value of the chosen item
        """
        return self.data.__getitem__(item)

    def __setitem__(self, item: int, value: float):
        """
        Sets the element at a certain position to a certain value.

        item: position
        value: new value
        """
        self.data.__setitem__(item, value)

    def __str__(self):
        """
        The string form of the matrix. Contains labels.

        returns: the prettified string
        """
        return tabulate(self.data, headers='keys', tablefmt='psql')

    def copy(self):
        """
        Copy current instance.

        returns: new instance
        """
        obj = type(self).__new__(self.__class__)
        obj.__dict__.update(self.__dict__)
        return obj

    def add(self, labels: List):
        """
        Add new labels and corresponding rows and columns.

        labels: labels to add
        """
        for label in labels:
            # Insert column.
            self.data[label] = self.init_value

            # Add label to the list of labels and increase size.
            self.labels.append(label)
            self.size = len(self.labels)

            # Insert row.
            self.data.loc[label] = [self.init_value] * (self.size)

    def drop(self, labels: List):
        """
        Drop labels and corresponding rows and columns.

        labels: labels to drop
        """
        self.data.drop(labels, inplace=True, axis=0)
        self.data.drop(labels, inplace=True, axis=1)
        for l in labels:
            self.labels.remove(l)
        self.size = len(self.labels)

    def argmax(self):
        """
        Returns the label belonging to the maximum value.
        
        returns the corresponding label
        """
        return (self.data.max(axis=0).idxmax(), self.data.max(axis=1).idxmax())

    def argmin(self, only_positive: bool = False):
        """
        Returns the label belonging to the minimum value.

        only_positive: look for minimum among the positive values
        returns the corresponding label
        """
        if only_positive:
            mask = (self.data <= 0)
            masked_data = self.data.mask(mask, np.inf, inplace=False)
            return (masked_data.min(axis=0).idxmin(), masked_data.min(axis=1).idxmin())
        return (self.data.min(axis=0).idxmin(), self.data.min(axis=1).idxmin())

    def posmax(self):
        """
        Returns the index belonging to the maximum value.

        returns the corresponding index
        """
        a, b = np.unravel_index(np.argmax(self.data.values, axis=None), self.data.values.shape)
        return b, a

    def posmin(self, only_positive: bool = False):
        """
        Returns the index belonging to the minimum value.

        only_positive: look for minimum among the positive values
        returns the corresponding index
        """
        if only_positive:
            mask = (self.data <= 0)
            masked_data = self.data.mask(mask, np.inf, inplace=False)
            a, b = np.unravel_index(np.argmin(masked_data.values, axis=None), masked_data.values.shape)
        else:
            a, b = np.unravel_index(np.argmin(self.data.values, axis=None), self.data.values.shape)
        return b, a
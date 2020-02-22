
from Bio.Phylo.BaseTree import Clade

class ParsimonyClade(Clade):
    
    def __init__(
        self,
        branch_length=None,
        name=None,
        clades=None,
        confidence=None,
        color=None,
        width=None,
        parent=None
    ):
        super(ParsimonyClade, self).__init__(branch_length, name, clades, confidence, color, width)
        self.parent = parent

    def __str__(self):
        return super(ParsimonyClade, self).__str__()

def remove_unprocessed(self, lst: list):
    """
    Remove clades that do not occur as many times as how many children they have.

    lst: list to process
    returns: new list
    """
    counts = Counter(lst)
    return [elem for elem in lst if counts[elem] == len(elem.clades)]

from Bio.Phylo import BaseTree
from Bio.Phylo import draw

parents = {}
clade0 = BaseTree.Clade(None, "Node0")
clade1 = BaseTree.Clade(None, "Node1")
clade2 = BaseTree.Clade(None, "Node2")
parents[clade1] = clade0
parents[clade2] = clade1
clade0.clades.append(clade1)
clade1.clades.append(clade2)
tree = BaseTree.Tree(clade0, rooted=False)
print(tree)
print(parents)

print()
parents[clade1].name = "N0"
print(tree)
print(parents)
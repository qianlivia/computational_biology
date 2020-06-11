from Bio.Phylo import BaseTree
from Bio.Phylo import draw

names = ["A", "B", "C", "D"]

clades = [BaseTree.Clade(None, name) for name in names]
clade1 = clades[0]
clade2 = clades[1]
print(clades)

inner_clade = BaseTree.Clade(None, "Inner0")
inner_clade.clades.append(clade1)
inner_clade.clades.append(clade2)

clades.append(inner_clade)
del clades[1]
del clades[0]
print(clades)

inner_clade = BaseTree.Clade(None, "Inner1")
inner_clade.clades.append(clades[0])
inner_clade.clades.append(clades[1])
clades.append(inner_clade)
del clades[1]
del clades[0]
print(clades)

inner_clade = BaseTree.Clade(None, "Inner2")
inner_clade.clades.append(clades[0])
inner_clade.clades.append(clades[1])
tree = BaseTree.Tree(inner_clade, rooted=True)

print(tree)

draw(tree)
from Bio.Phylo import BaseTree

names = ["A", "B", "C"]

clades = [BaseTree.Clade(None, name) for name in names]
clade1 = clades[0]
clade2 = clades[1]
print(clades)

inner_clade = BaseTree.Clade(None, "Inner0")
inner_clade.clades.append(clade1)
inner_clade.clades.append(clade2)

clades[0] = inner_clade
del clades[1]
print(clades)

inner_clade = BaseTree.Clade(None, "Inner1")
inner_clade.clades.append(clades[0])
inner_clade.clades.append(clades[1])
tree = BaseTree.Tree(inner_clade, rooted=False)

print(tree)
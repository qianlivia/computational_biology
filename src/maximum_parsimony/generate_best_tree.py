from Bio.Phylo import BaseTree, draw
from Bio.Phylo.BaseTree import Clade

# Coding

root = BaseTree.Clade(None, "")
clade_17735 = BaseTree.Clade(None, "")
clade_cat = BaseTree.Clade(None, "cat")
clade_8758 = BaseTree.Clade(None, "")
clade_7481 = BaseTree.Clade(None, "")
clade_4413 = BaseTree.Clade(None, "")
clade_3164 = BaseTree.Clade(None, "")
clade_2356 = BaseTree.Clade(None, "")
clade_horse = BaseTree.Clade(None, "horse")
clade_muntjak = BaseTree.Clade(None, "muntjak indian")
clade_2500 = BaseTree.Clade(None, "")
clade_gorilla = BaseTree.Clade(None, "gorilla")
clade_orangutan = BaseTree.Clade(None, "orangutan")
clade_owl_monkey = BaseTree.Clade(None, "owl monkey")
clade_lemur = BaseTree.Clade(None, "lemur")
clade_rat = BaseTree.Clade(None, "rat")
clade_mouse = BaseTree.Clade(None, "mouse")

root.clades.append(clade_17735)
root.clades.append(clade_cat)
clade_17735.clades.append(clade_8758)
clade_17735.clades.append(clade_7481)
clade_8758.clades.append(clade_4413)
clade_8758.clades.append(clade_2356)
clade_7481.clades.append(clade_horse)
clade_7481.clades.append(clade_muntjak)
clade_4413.clades.append(clade_3164)
clade_4413.clades.append(clade_lemur)
clade_2356.clades.append(clade_rat)
clade_2356.clades.append(clade_mouse)
clade_3164.clades.append(clade_2500)
clade_3164.clades.append(clade_owl_monkey)
clade_2500.clades.append(clade_gorilla)
clade_2500.clades.append(clade_orangutan)

clade_orangutan.branch_length = 2500
clade_gorilla.branch_length = 2500
clade_2500.branch_length = (3164 - 2500)
clade_owl_monkey.branch_length = 3164
clade_3164.branch_length = (4413 - 3164)
clade_lemur.branch_length = 4413
clade_4413.branch_length = (8758 - 4413)
clade_2356.branch_length = (8758 - 2356)
clade_rat.branch_length = 2356
clade_mouse.branch_length = 2356
clade_8758.branch_length = (17735 - 8758)
clade_7481.branch_length = (17735 - 7481)
clade_horse.branch_length = 7481
clade_muntjak.branch_length = 7481
clade_17735.branch_length = 18960 - 17735
clade_cat.branch_length = 18960

tree = BaseTree.Tree(root, rooted=True)

print(tree)

draw(tree)

##################################
# Non-coding

clade_131584 = BaseTree.Clade(None, "")
clade_120079 = BaseTree.Clade(None, "")
clade_93169 = BaseTree.Clade(None, "")
clade_55571 = BaseTree.Clade(None, "")
clade_35195 = BaseTree.Clade(None, "")
clade_22904 = BaseTree.Clade(None, "")
clade_28998 = BaseTree.Clade(None, "")
clade_horse = BaseTree.Clade(None, "horse")
clade_cow = BaseTree.Clade(None, "cow")
clade_gorilla = BaseTree.Clade(None, "gorilla")
clade_orangutan = BaseTree.Clade(None, "orangutan")
clade_owl_monkey = BaseTree.Clade(None, "owl monkey")
clade_lemur = BaseTree.Clade(None, "lemur")
clade_rat = BaseTree.Clade(None, "rat")
clade_mouse = BaseTree.Clade(None, "mouse")

clade_131584.clades.append(clade_120079)
clade_131584.clades.append(clade_mouse)
clade_120079.clades.append(clade_93169)
clade_120079.clades.append(clade_rat)
clade_93169.clades.append(clade_55571)
clade_93169.clades.append(clade_28998)
clade_55571.clades.append(clade_35195)
clade_55571.clades.append(clade_lemur)
clade_35195.clades.append(clade_22904)
clade_35195.clades.append(clade_owl_monkey)
clade_22904.clades.append(clade_gorilla)
clade_22904.clades.append(clade_orangutan)
clade_28998.clades.append(clade_cow)
clade_28998.clades.append(clade_horse)

clade_orangutan.branch_length = 22904
clade_gorilla.branch_length = 22904
clade_owl_monkey.branch_length = 35195
clade_lemur.branch_length = 55571
clade_cow.branch_length = 28998
clade_horse.branch_length = 28998
clade_rat.branch_length = 120079
clade_mouse.branch_length = 131584
clade_22904.branch_length = 35195 - 22904
clade_35195.branch_length = 55571 - 35195
clade_55571.branch_length = 93169 - 55571
clade_93169.branch_length = 120079 - 93169
clade_120079.branch_length = 131584 - 120079
clade_28998.branch_length = 93169 - 28998

tree = BaseTree.Tree(clade_131584, rooted=True)

print(tree)

draw(tree)

##################################
# Non-coding 2

clade_200477 = BaseTree.Clade(None, "")
clade_186005 = BaseTree.Clade(None, "")
clade_137114 = BaseTree.Clade(None, "")
clade_100069 = BaseTree.Clade(None, "")
clade_62477 = BaseTree.Clade(None, "")
clade_27787 = BaseTree.Clade(None, "")
clade_17990 = BaseTree.Clade(None, "")
clade_tetra = BaseTree.Clade(None, "tetra")
clade_shrew = BaseTree.Clade(None, "shrew")
clade_gorilla = BaseTree.Clade(None, "gorilla")
clade_chimp = BaseTree.Clade(None, "chimp")
clade_squirrel_monkey = BaseTree.Clade(None, "squirrel monkey")
clade_guinea_pig = BaseTree.Clade(None, "guinea pig")
clade_baboon = BaseTree.Clade(None, "baboon")
clade_galago = BaseTree.Clade(None, "galago")

clade_200477.clades.append(clade_186005)
clade_200477.clades.append(clade_galago)
clade_186005.clades.append(clade_137114)
clade_186005.clades.append(clade_27787)
clade_137114.clades.append(clade_100069)
clade_137114.clades.append(clade_guinea_pig)
clade_100069.clades.append(clade_62477)
clade_100069.clades.append(clade_squirrel_monkey)
clade_62477.clades.append(clade_tetra)
clade_62477.clades.append(clade_shrew)
clade_27787.clades.append(clade_baboon)
clade_27787.clades.append(clade_17990)
clade_17990.clades.append(clade_gorilla)
clade_17990.clades.append(clade_chimp)

clade_tetra.branch_length = 62477
clade_shrew.branch_length = 62477
clade_squirrel_monkey.branch_length = 100069
clade_62477.branch_length = 100069 - 62477
clade_100069.branch_length = 137114 - 100069
clade_guinea_pig.branch_length = 137114
clade_137114.branch_length = 186005 - 137114
clade_186005.branch_length = 200477 - 186005
clade_27787.branch_length = 186005 - 27787
clade_baboon.branch_length = 27787
clade_gorilla.branch_length = 17990
clade_chimp.branch_length = 17990
clade_galago.branch_length = 200477
clade_17990.branch_length = 27787 - 17990

tree = BaseTree.Tree(clade_200477, rooted=True)

print(tree)

draw(tree)
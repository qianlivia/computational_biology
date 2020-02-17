
from utils import Matrix, Vector
import numpy as np
import pandas as pd

mat = Matrix(["A", "B", "C"])
mat["A"]["A"] = 1
mat["A"]["B"] = -13
mat["C"]["C"] = -543
mat.add("D")
mat["D"]["C"] = -1432
print(mat)
print(mat.argmin())
print(mat.posmin())
mat.drop(["A"])
print(mat)
mat2 = mat.copy()
print("Copy")
print(mat2)
print(mat.labels, mat2.labels, mat.size, mat2.size)


vec = Vector(["A", "B", "C"])
vec["A"] = 1
vec["B"] = -13
print(vec)
print(vec.argmin())
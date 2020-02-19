
from tools.utils import Matrix, Vector
import numpy as np
import pandas as pd

mat = Matrix(["A", "B", "C"], init_value=1)
mat.add("D")
mat["A"]["B"] = 13
mat["A"]["C"] = 1
mat["A"]["D"] = 10
mat["B"]["C"] = 543
mat["B"]["D"] = 1
mat["C"]["D"] = 1432
print(mat)
print(mat.argmin(only_positive=True))
print(mat.posmin(only_positive=True))
mat.drop(["A"])
print(mat)
mat2 = mat.copy()
print("Copy")
print(mat2)
print(mat.labels, mat2.labels, mat.size, mat2.size)


vec = Vector(["A", "B", "C"], init_value=np.inf)
vec["A"] = 1
vec["B"] = -13
vec["D"] = 30
print(vec)
print(vec.argmin())
vec.add("D")
vec.drop(["B"])
print(vec)

from utils import Matrix, Vector
import numpy as np
import pandas as pd

mat = Matrix(["A", "B", "C"])
mat["A"]["A"] = 1
mat["A"]["B"] = -13
print(mat)
print(mat.argmin())
print(mat.posmin())
mat.drop_columns(["A"])
mat.drop_rows(["A"])
print(mat)

vec = Vector(["A", "B", "C"])
vec["A"] = 1
vec["B"] = -13
print(vec)
print(vec.argmin())

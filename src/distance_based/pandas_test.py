
from distance_calculator import Matrix, Vector
import numpy as np

mat = Matrix(["A", "B", "C"])
mat["A"]["A"] = 1
mat["A"]["B"] = -13
print(mat)
print(mat.argmin())
print(mat.posmin())

vec = Vector(["A", "B", "C"])
vec["A"] = 1
vec["B"] = -13
print(vec)
print(vec.argmin())
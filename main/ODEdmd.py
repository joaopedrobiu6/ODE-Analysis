import matplotlib.pyplot as plt
import numpy as np

from pydmd import DMD

file_in = open("data1.csv", "r")
r = file_in.readlines()
file_in.close()
#l = 0
#x = []
#while l < len(r):
#    floating = float(r[l])
#    x.append(floating)
#    l += 1

print(r)

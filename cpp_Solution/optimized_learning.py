import numpy as np
import time
import random
from tqdm import tqdm

# import sys
# sys.path.insert(0, "C:\\Users\\samia\\a_doMath\\domath2021\\cpp_Solution")
import pysol

sol_obj = pysol.PySolution()

#First Parameter is number of iterations per starting point
#Second Parameter is number of starting points, if starting point is not under threshold, then point is not optimized
#Third Parameter is target value that a starting point needs to be less than for optimization to begin

c, params = sol_obj.optimize(100, 1000, 100)

print("Best c Value ", c)
print("Parameters ", params)
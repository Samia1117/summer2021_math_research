cdef extern from "Solution.cpp":
	pass

# Declare the class with cdef
cdef extern from "Solution.h":
	cdef struct opVal:
		float c
		float** params
	cdef cppclass Solution:
		Solution() except +
		float** create_random_params(float** centers, float width)
		float i_fast(float** params, int intervals)
		opVal optimize(int iterations, int num_starts, float target)
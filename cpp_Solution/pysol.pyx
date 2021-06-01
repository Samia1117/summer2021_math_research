# distutils: language = c++

import numpy as np
cimport numpy as np

np.import_array()

from Solution cimport Solution
from Solution cimport opVal
from array_pointer cimport *

# Create a Cython extension type which holds a C++ instance
# as an attribute and create a bunch of forwarding methods
# Python extension type.

cdef class PySolution:
	cdef Solution sol  # Hold a C++ instance which we're wrapping

	def __cinit__(self):
		self.sol = Solution()

	cpdef create_random_params(self, centers, float width):
		cdef float** inputCenters = arr_to_ptr_ptr(2)

		cdef float* centers_a = arr_to_ptr(5)
		cdef float* centers_b = arr_to_ptr(5)
		cdef int k
		for k in range(5):
			centers_a[k] = centers[0][k]
			centers_b[k] = centers[1][k]
		inputCenters[0] = centers_a
		inputCenters[1] = centers_b

		cdef float** arr = self.sol.create_random_params(inputCenters, width)
		cdef np.ndarray[float, ndim = 2] ret
		ret = np.zeros((2,5)).astype(np.single)

		cdef int i,j
		for i in range(2):
			for j in range(5):
				ret[i][j] = arr[i][j]	
		return ret
	
	cpdef float i_fast(self, params, int intervals):
		cdef float** ret = arr_to_ptr_ptr(2)

		cdef float* p_a = arr_to_ptr(5)
		cdef float* p_b = arr_to_ptr(5)
		cdef int i,j
		for i in range(6):
			p_a[i] = params[0][i]
			p_b[i] = params[1][i]
		ret[0] = p_a
		ret[1] = p_b
		
		return self.sol.i_fast(ret, intervals)
	
	cpdef optimize(self, int iterations, int num_starts, float target):
		cdef opVal ret = self.sol.optimize(iterations, num_starts, target)
		
		cdef np.ndarray[float, ndim = 2] params
		params = np.zeros((2,5)).astype(np.single)

		cdef int i,j
		for i in range(2):
			for j in range(5):
				params[i][j] = ret.params[i][j]	

		return ret.c, params
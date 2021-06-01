cdef extern from "array_pointer.cpp":
	pass

#needed functions for pysol.pyx file
cdef extern from "array_pointer.h":
	float* arr_to_ptr(int m)
	float** arr_to_ptr_ptr(int m)
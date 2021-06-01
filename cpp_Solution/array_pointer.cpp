#include "array_pointer.h"

float* arr_to_ptr(int m) {
	float* ret = new float[m];
	return ret;
}

float** arr_to_ptr_ptr(int m) {
	float** ret = new float*[m];
	return ret;
}
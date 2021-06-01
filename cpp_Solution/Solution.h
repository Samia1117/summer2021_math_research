#ifndef SOLUTION_H
#define SOLUTION_H

struct opVal;

class Solution {
	public:
		Solution();
		float curWidth;
		float** curCenters;
		void setWidth(float width);
		void setCenters(float** centers);
		float** create_random_params(float** centers, float width);
		float eval_f(float* x_params, float x);
		float eval_f_prime(float* x_params, float x);
		float i_fast(float** params, int intervals);
		opVal optimize(int iterations, int num_starts, float target);
};

#endif
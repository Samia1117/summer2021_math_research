#include "Solution.h"
#include <iostream>

#include <ctime>
#include <cstdlib>
#include <math.h>

struct opVal {
	float c;
	float** params;
};

Solution::Solution() {
	srand (static_cast <unsigned> (time(0)));
	float** defc = new float*[1];
	for (int i = 0; i<2; i++) {
		defc[i] = new float[12]{0,0,0,0,0,0,0,0,0,0,0,0};
	}
	this->curWidth = 1.;
	this->curCenters = defc;
}

void Solution::setWidth(float width) {
	this->curWidth = width;
}

void Solution::setCenters(float** centers) {
	this->curCenters = centers;
}

float** Solution::create_random_params(float** centers, float width) {
	float** ret = new float*[2];
	for (int i = 0; i < 2; i++) {
		ret[i] = new float[12];
		for (int j = 0; j < 12; j++) {
			float LO = centers[i][j] - width;
			float HI = centers[i][j] + width;
			ret[i][j] = LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
		}
	}
	return ret;
};
// f(x) = a_0*cos(x) + a_1*sin(x) + a_2*cos(2*x) + a_3*sin(2*x) + a_4*cos(3*x) 
// + a_5*sin(3*x) + a_6*cos(4*x) + a_7*sin(4*x) + a_8*cos(5*x) + a_9*sin(5*x)
float Solution::eval_f(float* x_params, float x){
	float f = x_params[0]*cos(x) + x_params[1]*sin(x) + x_params[2]*cos(2*x) + x_params[3]*sin(2*x) +
	x_params[4]*cos(3*x) + x_params[5]*sin(3*x) + x_params[6]*cos(4*x) + x_params[7]*sin(4*x) 
	+ x_params[8]*cos(5*x) + x_params[9]*sin(5*x) + (x_params[10]*cos(6*x)) + x_params[11]*sin(6*x);
    return f;
};
// f'(x) = -a_0*sin(x) + a_1*cos(x) + -2*a_2*sin(2*x) + 2*a_3*cos(2*x) + -3*a_4*sin(3*x) 
// + 3*a_5*cos(3*x) + -4*a_6*sin(4*x) + 4*a_7*cos(4*x) + -5*a_8*sin(5*x) + 5*a_9*sin(5*x)
float Solution::eval_f_prime(float* x_params, float x){
    float f_prime = -x_params[0]*sin(x) + x_params[1]*cos(x) + (-2*x_params[2]*sin(2*x)) + 2*x_params[3]*cos(2*x) +
	(-3*x_params[4]*sin(3*x)) + 3*x_params[5]*cos(3*x) + (-4*x_params[6]*sin(4*x)) + 4*x_params[7]*cos(4*x) +
	(-5*x_params[8]*sin(5*x)) + 5*x_params[9]*cos(5*x) + (-6*x_params[10]*sin(6*x)) + 6*x_params[11]*cos(6*x);
	return f_prime;
};

float Solution::i_fast(float** params, int intervals) {
	float* x_params = params[0];
	float PI = 3.14159;
    float* vals = new float[intervals];
	float x;

    float rhs = 0;
	float lhs = 0;

    for (int i=0; i<intervals; i++){
        vals[i] = -PI + (2*PI*i)/(intervals-1);
    }
    for (int i = 0; i<intervals; i++) {
        float x = vals[i];
		float f_prime = eval_f_prime(x_params, x);
		float f = eval_f(x_params, x);

		rhs += pow(f_prime, 3); // det_u*del_u^2
        lhs += pow((pow(f_prime, 2) +1), 2) + ((2/3)*pow(f, 2)*(pow(f_prime, 2)+1)) + pow(f, 4); //del_u^4
    }
    if (rhs == 0) {
		return 0;
	}
	return lhs/rhs;
};

opVal Solution::optimize(int iterations, int num_starts, float target) {
	std::cout << "Calling Optimize Function" << std::endl;

	float MAX_C = 9999999;
	float** maxvals = new float*[2];
	for (int i = 0; i<2; i++) {
		maxvals[i] = new float[12]{0,0,0,0,0,0,0,0,0,0,0,0};
	}

	for (int i = 0; i<num_starts; i++) {
		if (i % 500 == 0) {
			std::cout << "Finished " << i << " out of " << num_starts << " points." << std::endl;
		}
		float** resetCenters = new float*[2];
		for (int i = 0; i<2; i++) {
			resetCenters[i] = new float[12]{0,0,0,0,0,0,0,0,0,0,0,0};
		}
		setCenters(resetCenters);
		setWidth(1);

		float** temp_max_params = create_random_params(this->curCenters, this->curWidth);
		float test_c = i_fast(temp_max_params,500);
		float temp_max_c = 0;
		if (test_c > 0 && test_c <= target + 2) {
			temp_max_c = i_fast(temp_max_params,1500);
		}
		else {
			continue;
		}
		//std::cout << "First c Value " << temp_max_c << std::endl;
		if (temp_max_c > 0 && temp_max_c <= target) {
			std::cout << "Beginning to Optimize Point: " << temp_max_c << std::endl;
			setCenters(temp_max_params);
		
			for (int j = 0; j<iterations; j++) {
				float** params = create_random_params(this->curCenters, this->curWidth);
				float temp_c = i_fast(params, 1000);
				if (temp_c <= 0)	continue;
				else {
					if (temp_c < temp_max_c) {
						float maxdiff = 0;
						for (int k = 0; k < 2; k++) {
							for (int l = 0; l<12; l++) {
								if(abs(temp_max_params[k][l] - params[k][l] > maxdiff))	maxdiff = abs(temp_max_params[k][l] - params[k][l] > maxdiff);
							}
						}
						temp_max_c = temp_c;
						temp_max_params = params;
						setCenters(params);
						setWidth(maxdiff);
					}
				}
			}
			std::cout << "Optimized c value " << temp_max_c << std::endl;

			if (temp_max_c < MAX_C) {
				MAX_C = temp_max_c;
				maxvals = temp_max_params;
			}
		}
	}
	opVal ret;
	ret.c = MAX_C;
	ret.params = maxvals;
	return ret;
}




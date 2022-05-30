/*	TestRegression.cpp
	Ahmed M. Hussein (amhussein4@gmail.com)
	05/05/2022

Copyright (c) 2013 Ahmed M. Hussein (amhussein4@gmail.com)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#include "stdio.h"
#include "stdlib.h"

extern "C"
{
	#include "BLAS.h"
	#include "Random.h"
	#include "Regression.h"
}

void TestRegression(unsigned int data_size)
{
	// sample data from a noisy linear model
	unsigned int input_size = 8;
	unsigned int output_size = 6;
	Matrix* transformation = create_matrix(output_size,input_size);
	Matrix* bias = create_matrix(output_size,1);
	mat_rand(transformation,-10.0,10.0,0);
	mat_rand(bias,-4.0,4.0,0);
	Matrix* input = create_matrix(input_size,1);
	Matrix* output = create_matrix(output_size,1);
	Matrix* noise = create_matrix(output_size,1);
	RegressionData* data = create_regdata(data_size,input_size,output_size);
	for(unsigned int i = 0 ; i < data_size ; i++)
	{
		mat_rand(input,-100.0,100.0,0);
		mat_nrand(noise,0.0,3.0);
		mat_mult(transformation,input,output);
		mat_inc(output,bias);
		mat_inc(output,noise);
		reg_set_entry(data,i,input->entries,output->entries);
	}
	RegressionModel* model = create_regmodel(input_size,output_size);
	reg_fit(data,model);
	//printf("orig weights\n");
	//mat_print(model->weights);
	Matrix* weights_t = create_matrix(1,1);
	mat_transpose(model->weights,weights_t);
	Matrix* affine_t = create_matrix(output_size,input_size + 1);
	for(unsigned int i = 0 ; i < output_size ; i++)
	{
		mat_set(affine_t,i,0,mat_get(bias,i,0));
		dcopy(input_size,1,transformation->entries + i*input_size,1,affine_t->entries + i*(input_size + 1) + 1);
	}
	Matrix* diff = create_matrix(1,1);
	mat_diff(weights_t,affine_t,diff);
	//mat_print(diff);
	printf("error norm2 = %e\n",mat_norm2(diff));
	double* test_inputs = new double[input_size];
	double* test_outputs = new double[output_size];

	for(unsigned int i = 0 ; i < input_size ; i++)
	{
		test_inputs[i] = rand_uniform_interval(-100.0,100.0);
		input->entries[i] = test_inputs[i];
	}
	mat_mult(transformation,input,output);
	mat_inc(output,bias);
	mat_print(output);
	//mat_print(transformation);
	//mat_print(bias);
	reg_estimate(model,test_inputs,test_outputs);
	for(unsigned int i = 0 ; i < output_size ; i++)
	{
		printf("%u : %e\n",i,test_outputs[i]);
	}
	delete [] test_inputs;
	delete [] test_outputs;
	destroy_matrix(transformation);
	destroy_matrix(bias);
	destroy_matrix(input);
	destroy_matrix(output);
	destroy_matrix(noise);
	destroy_matrix(weights_t);
	destroy_matrix(affine_t);
	destroy_matrix(diff);
	destroy_regdata(data);
	destroy_regmodel(model);
}

int main(int argc,char** argv)
{
	if(argc < 2)
	{
		printf("error: missing run arguments\n");
		printf("usage:testregression data_size\n");
		return 1;
	}
	unsigned int data_size = atoi(argv[1]);
	TestRegression(data_size);
	return 0;
}


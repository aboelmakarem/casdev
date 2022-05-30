/*	Regression.c
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

#include "Regression.h"
#include "BLAS.h"
#include "stdlib.h"

RegressionData* create_regdata(unsigned int data_size,unsigned int input_size,unsigned int output_size)
{
	RegressionData* data = (RegressionData*)malloc(sizeof(RegressionData));
	data->inputs = create_matrix(data_size,input_size);
	data->outputs = create_matrix(data_size,output_size);
	return data;
}
void reset_regdata(RegressionData* data)
{
	if(data == 0)					return;
	if(data->inputs != 0)			destroy_matrix(data->inputs);
	if(data->outputs != 0)			destroy_matrix(data->outputs);
	data->inputs = 0;
	data->outputs = 0;
}
void destroy_regdata(RegressionData* data)
{
	if(data == 0)					return;
	reset_regdata(data);
	free(data);
}
unsigned int reg_data_size(const RegressionData* data){return data->inputs->rows;}
unsigned int reg_input_size(const RegressionData* data){return data->inputs->columns;}
unsigned int reg_output_size(const RegressionData* data){return data->outputs->columns;}
int reg_set_entry(RegressionData* data,unsigned int index,const double* inputs,const double* outputs)
{
	if(index >= data->inputs->rows)			return 0;
	dcopy(data->inputs->columns,1,inputs,1,data->inputs->entries + index*data->inputs->columns);
	dcopy(data->outputs->columns,1,outputs,1,data->outputs->entries + index*data->outputs->columns);
	return 1;
}
int reg_set_input(RegressionData* data,const Matrix* inputs)
{
	if(inputs->rows != data->inputs->rows)					return 0;
	if(inputs->columns != data->inputs->columns)			return 0;
	copy_matrix(inputs,data->inputs);
	return 1;
}
int reg_set_output(RegressionData* data,const Matrix* outputs)
{
	if(outputs->rows != data->outputs->rows)				return 0;
	if(outputs->columns != data->outputs->columns)			return 0;
	copy_matrix(outputs,data->outputs);
	return 1;
}
int reg_set(RegressionData* data,const Matrix* inputs,const Matrix* outputs)
{
	if(!reg_set_input(data,inputs))			return 0;
	if(!reg_set_output(data,outputs))		return 0;
	return 1;
}

RegressionModel* create_regmodel(unsigned int input_size,unsigned int output_size)
{
	RegressionModel* model = (RegressionModel*)malloc(sizeof(RegressionModel));
	model->weights = create_matrix(input_size + 1,output_size);
	return model;
}
void reset_regmodel(RegressionModel* model)
{
	if(model == 0)					return;
	if(model->weights != 0)			destroy_matrix(model->weights);
	model->weights = 0;
}
void destroy_regmodel(RegressionModel* model)
{
	if(model == 0)					return;
	reset_regmodel(model);
	free(model);
}
void reg_estimate(const RegressionModel* model,const double* inputs,double* outputs)
{
	// subtract 1 to account for the bias term in the weight matrix
	unsigned int input_size = model->weights->rows - 1;
	unsigned int output_size = model->weights->columns;
	for(unsigned int i = 0 ; i < output_size ; i++)
	{
		outputs[i] = model->weights->entries[i] + ddot(input_size,1,inputs,output_size,model->weights->entries + i + output_size);
	}
}
int reg_fit(const RegressionData* data,RegressionModel* model)
{
	unsigned int data_size = data->outputs->rows;
	if(data->inputs->rows != data_size)				return 0;
	unsigned int input_size = data->inputs->columns;
	unsigned int output_size = data->outputs->columns;
	// weight_count is the number of weights per output dimension
	// an extra weight is added for the bias term, it is always multiplied by 1
	unsigned int weight_count = input_size + 1;
	// create phi matrix which has 1 extra input of 1 for bias term
	// all other entries of the phi matrix are similar to the input matrix
	Matrix* phi = create_matrix(data_size,weight_count);
	for(unsigned int i = 0 ; i < data_size ; i++)
	{
		mat_set(phi,i,0,1.0);
		dcopy(input_size,1,data->inputs->entries + i*input_size,1,phi->entries + i*weight_count + 1);
	}
	// build and solve the normal equation: phi_t phi = phi^t t where phi is the matrix of inputs 
	// (an n x m matrix where n is the number of data points and m is the number of inputs)
	// t is the output matrix (an n x d matrix) where d is the output vector dimensions
	Matrix* phi_t = create_matrix(weight_count,data_size);
	mat_transpose(phi,phi_t);
	// normal equation system matrix and right hand side
	Matrix* A = create_matrix(weight_count,weight_count);
	mat_mult(phi_t,phi,A);
	Matrix* rhs = create_matrix(weight_count,output_size);
	mat_mult(phi_t,data->outputs,rhs);
	//mat_solve_ge(A,rhs,model->weights);
	// solution using GMRES requires passing each right hand side vector in a separate call
	Matrix* b = create_matrix(weight_count,1);
	Matrix* x = create_matrix(weight_count,1);
	for(unsigned int i = 0 ; i < output_size ; i++)
	{
		dcopy(weight_count,output_size,rhs->entries + i,1,b->entries);
		mat_solve_gmres(A,b,x);
		dcopy(weight_count,1,x->entries,output_size,model->weights->entries + i);
	}
	destroy_matrix(phi);
	destroy_matrix(phi_t);
	destroy_matrix(A);
	destroy_matrix(rhs);
	destroy_matrix(b);
	destroy_matrix(x);
	return 1;
}


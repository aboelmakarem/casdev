/*	Regression.h
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

#ifndef REGRESSION_H_
#define REGRESSION_H_

#include "Matrix.h"

typedef struct RegressionData
{
	Matrix* inputs;
	Matrix* outputs;
} RegressionData;

typedef struct RegressionModel
{
	Matrix* weights;
} RegressionModel;

// Regression data
RegressionData* create_regdata(unsigned int data_size,unsigned int input_size,unsigned int output_size);
void reset_regdata(RegressionData* data);
void destroy_regdata(RegressionData* data);
unsigned int reg_data_size(const RegressionData* data);
unsigned int reg_input_size(const RegressionData* data);
unsigned int reg_output_size(const RegressionData* data);
int reg_set_entry(RegressionData* data,unsigned int index,const double* inputs,const double* outputs);
int reg_set_input(RegressionData* data,const Matrix* inputs);
int reg_set_output(RegressionData* data,const Matrix* outputs);
int reg_set(RegressionData* data,const Matrix* inputs,const Matrix* outputs);
// Regression model
RegressionModel* create_regmodel(unsigned int input_size,unsigned int output_size);
void reset_regmodel(RegressionModel* model);
void destroy_regmodel(RegressionModel* model);
void reg_estimate(const RegressionModel* model,const double* inputs,double* outputs);
// Model fitting
int reg_fit(const RegressionData* data,RegressionModel* model);


#endif


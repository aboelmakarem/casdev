/*	BLAS.h
	Ahmed M. Hussein (amhussein4@gmail.com)
	04/27/2022

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

#ifndef BLAS_H_
#define BLAS_H_

void dcopy(unsigned int size,unsigned int x_inc,const double* x,unsigned int y_inc,double* y);
void dswap(unsigned int size,unsigned int x_inc,double* x,unsigned int y_inc,double* y);
void dscal(unsigned int size,unsigned int x_inc,double* x,double factor);
void daxpy(unsigned int size,unsigned int x_inc,const double* x,double factor,unsigned int y_inc,double* y);
void drotg(double a,double b,double* angle_cosine,double* angle_sine);
void drot(unsigned int size,unsigned int x_inc,double* x,unsigned int y_inc,double* y,double angle_cosine,double angle_sine);
double dasum(unsigned int size,unsigned int x_inc,const double* x);
double ddot(unsigned int size,unsigned int x_inc,const double* x,unsigned int y_inc,const double* y);
double dnrm2(unsigned int size,unsigned int x_inc,const double* x);
double dnrm(unsigned int size,unsigned int x_inc,const double* x);
unsigned int idamax(unsigned int size,unsigned int x_inc,const double* x);
void dgemv(unsigned int row_count,unsigned int column_count,double alpha,const double* matrix,unsigned int x_inc,const double* x,double beta,unsigned int y_inc,double* y,int operation);
void dtrmv(unsigned int size,const double* matrix,unsigned int x_inc,const double* x,unsigned int y_inc,double* y,int matrix_type);
void dtrsv(unsigned int size,const double* matrix,unsigned int x_inc,double* x,int matrix_type);
void ddgmv(unsigned int size,const double* matrix,unsigned int x_inc,const double* x,unsigned int y_inc,double* y,int matrix_type);
int test_BLAS(unsigned int size,int comprehensive,double tolerance,unsigned int paranoia);
// testing functions
int test_dcopy(unsigned int test_size,unsigned int test_paranoia,double tolerance);
int test_dswap(unsigned int test_size,unsigned int test_paranoia,double tolerance);
int test_dscal(unsigned int test_size,unsigned int test_paranoia,double tolerance);
int test_daxpy(unsigned int test_size,unsigned int test_paranoia,double tolerance);
int test_drotg(unsigned int test_paranoia,double tolerance);
int test_drot(unsigned int test_paranoia,double tolerance);
int test_dasum();
int test_ddot(unsigned int test_size,unsigned int test_paranoia,double tolerance);
int test_dnrm2(unsigned int test_size,unsigned int test_paranoia,double tolerance);
int test_dnrm(unsigned int test_size,unsigned int test_paranoia,double tolerance);
int test_idamax();
int test_dgemv(unsigned int test_size,unsigned int test_paranoia,double tolerance);
int test_dtrmv(unsigned int test_size,unsigned int test_paranoia,double tolerance);
int test_dtrsv(unsigned int test_size,unsigned int test_paranoia,double tolerance);
int test_ddgmv(unsigned int test_size,unsigned int test_paranoia,double tolerance);
void vector_random_populate(unsigned int size,double* x);
void matrix_random_populate(unsigned int row_count,unsigned int column_count,double* matrix,int non_singular);
void make_upper_triangular(unsigned int size,double* matrix);
void make_lower_triangular(unsigned int size,double* matrix);
void make_unit_diagonal(unsigned int size,double* matrix);
void set_value(unsigned int size,double* x,double value);
int compare_vectors(unsigned int size,unsigned int x_inc,const double* x,unsigned int y_inc,const double* y,double y_factor,unsigned int y_shift_inc,const double* y_shift,double tolerance);

#endif


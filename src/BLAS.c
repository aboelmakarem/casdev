/*	BLAS.c
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

#include "BLAS.h"
#include "math.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"

// BLAS parameters
double BLAS_tolerance = 1.0e-9;
// BLAS functions
void dcopy(unsigned int size,unsigned int x_inc,const double* x,unsigned int y_inc,double* y)
{
	// This function copies N elements of the array X to the array Y. The increments x_inc and y_inc 
	// are used so that the copied elements can be spaced at different intervals. This enables copying 
	// matrix rows into columns and the like regardless whether the matrix stores its data 
	// row-wise or column wise. 
	if(size == 0)				return;
	if((x_inc == 1) && (y_inc == 1))
	{
		// if we are doing a continuous, non strided, copy, do a fast block memory copy
		memcpy(y,x,size*sizeof(double));
	}
	else
	{
		unsigned int x_index = 0;
		unsigned int y_index = 0;
		for(unsigned int i = 0 ; i < size ; i++)
		{
			y[y_index] = x[x_index];
			x_index += x_inc;
			y_index += y_inc;
		}
	}
}
void dswap(unsigned int size,unsigned int x_inc,double* x,unsigned int y_inc,double* y)
{
	// This function swaps N elements of the array X and Y arrays. The increments x_inc and y_inc 
	// are used so that the swapped elements can be spaced at different intervals. This enables swapping 
	// matrix rows into columns and the like regardless whether the matrix stores its data 
	// row-wise or column wise. 
	if(size == 0)				return;
	double temp = 0.0;
	if((x_inc == 1) && (y_inc == 1))
	{
		// unroll the loops 3 at a time
		unsigned int step = 3;
		unsigned int i = 0;
		unsigned int m = size%step;
		for(i = 0 ; i < m ; i++)
		{
			temp = x[i];
			x[i] = y[i];
			y[i] = temp;
		}
		i = m;
		while(i < size)
		{
			temp = x[i];
			x[i] = y[i];
			y[i] = temp;
			i++;
			temp = x[i];
			x[i] = y[i];
			y[i] = temp;
			i++;
			temp = x[i];
			x[i] = y[i];
			y[i] = temp;
			i++;
		}
	}
	else
	{
		unsigned int x_index = 0;
		unsigned int y_index = 0;
		for(unsigned int i = 0 ; i < size ; i++)
		{
			temp = x[x_index];
			x[x_index] = y[y_index];
			y[y_index] = temp;
			x_index += x_inc;
			y_index += y_inc;
		}
	}
}
void dscal(unsigned int size,unsigned int x_inc,double* x,double factor)
{
	// This function scales N elements of the array X spaced x_inc apart by the factor Factor.
	if(size == 0)								return;
	if(x_inc == 0)								return;
	if(fabs(factor - 1.0) < BLAS_tolerance)		return;
	if(x_inc == 1)
	{
		// unroll the loops 5 at a time
		unsigned int step = 5;
		unsigned int i = 0;
		unsigned int m = size%step;
		for(i = 0 ; i < m ; i++)
		{
			x[i] = x[i]*factor;
		}
		i = m;
		while(i < size)
		{
			x[i] = x[i]*factor;
			x[i + 1] = x[i + 1]*factor;
			x[i + 2] = x[i + 2]*factor;
			x[i + 3] = x[i + 3]*factor;
			x[i + 4] = x[i + 4]*factor;
			i += step;
		}
	}
	else
	{
		unsigned int x_index = 0;
		for(unsigned int i = 0 ; i < size ; i++)
		{
			x[x_index] = x[x_index]*factor;
			x_index += x_inc;
		}
	}
}
void daxpy(unsigned int size,unsigned int x_inc,const double* x,double factor,unsigned int y_inc,double* y)
{
	// This function scales N elements of the array x spaced x_inc apart by the given factor then adds 
	// the scaled array to another array y with elements spaced y_inc apart and stores the result in Y.
	if(size == 0)								return;
	if(fabs(factor) < BLAS_tolerance)			return;
	if((x_inc == 1) && (y_inc == 1))
	{
		// unroll the loops 4 at a time
		unsigned int step = 4;
		unsigned int i = 0;
		unsigned int m = size%step;
		for(i = 0 ; i < m ; i++)
		{
			y[i] += x[i]*factor;
		}
		i = m;
		while(i < size)
		{
			y[i] += x[i]*factor;
			y[i + 1] += x[i + 1]*factor;
			y[i + 2] += x[i + 2]*factor;
			y[i + 3] += x[i + 3]*factor;
			i += step;
		}
	}
	else
	{
		unsigned int x_index = 0;
		unsigned int y_index = 0;
		for(unsigned int i = 0 ; i < size ; i++)
		{
			y[y_index] += x[x_index]*factor;
			x_index += x_inc;
			y_index += y_inc;
		}
	}
}
void drotg(double a,double b,double* angle_cosine,double* angle_sine)
{
	// This function computes the Givens rotation that would transform the 
	// vector [A B]^T to the vector [R 0]^T. This can be done by premultiplying 
	// the vector by G where G is an 2x2 orthonormal matrix of the form [c s ; -s c]
	// and its effect is to rotate the vector by some angle theta in its plane.
	// The cosine and sine of the rotation angle are stored in angle_cosine and angle_sine respectively.
	// The solution for C and S are C = A/R and S = B/R where R = sqrt(A^2 + B^2). 
	if(fabs(b) < BLAS_tolerance)
	{
		// if B is already zero, no rotation is needed
		(*angle_cosine) = 1.0;
		(*angle_sine) = 0.0;
		return;
	}
	// for stability, divide by the larger of the two vector entries
	double t = 0.0;
	double r = 0.0;
	if(fabs(b) > fabs(a))
	{
		t = a/b;
		r = 1.0/sqrt(1.0 + t*t);
		(*angle_sine) = r;
		(*angle_cosine) = r*t;
	}
	else
	{
		t = b/a;
		r = 1.0/sqrt(1.0 + t*t);
		(*angle_cosine) = r;
		(*angle_sine) = r*t;
	}
}
void drot(unsigned int size,unsigned int x_inc,double* x,unsigned int y_inc,double* y,double angle_cosine,double angle_sine)
{
	// This function applies a Givens rotation described by the cosine and the sine of an angle 
	// (C and S) on the two arrays X and Y both of size N with elements spaced x_inc and y_inc apart. 
	// Since applying Givens rotation to a matrix amounts to modifying two of its rows, a typical 
	// use case for this function is to pass the arrays of the two rows that will be modified
	// and they will be modified in place by this function. This is why the indices of the two 
	// rows to be modified is immaterial because the actual arrays, not their indices, will be passed.
	// The specification of x_inc and y_inc allows passing any rows directly even if the matrix is 
	// stored column-wise. 
	if(size == 0)				return;
	double temp = 0.0;
	if((x_inc == 1) && (y_inc == 1))
	{
		for(unsigned int i = 0 ; i < size ; i++)
		{
			temp = angle_cosine*x[i] + angle_sine*y[i];
			y[i] = angle_cosine*y[i] - angle_sine*x[i];
			x[i] = temp;
		}
	}
	else
	{
		unsigned int x_index = 0;
		unsigned int y_index = 0;
		for(unsigned int i = 0 ; i < size ; i++)
		{
			temp = angle_cosine*x[x_index] + angle_sine*y[y_index];
			y[y_index] = angle_cosine*y[y_index] - angle_sine*x[x_index];
			x[x_index] = temp;
			x_index += x_inc;
			y_index += y_inc;
		}
	}
}
double dasum(unsigned int size,unsigned int x_inc,const double* x)
{
	// This function computes the sum of N absolute values of the array X entries spaced x_inc 
	// apart.
	double sum = 0.0;
	if(size == 0)				return sum;
	if(x_inc == 0)				return sum;
	if(x_inc == 1)
	{
		// unroll the loops 6 at a time
		unsigned int step = 6;
		unsigned int i = 0;
		unsigned int m = size%step;
		for(i = 0 ; i < m ; i++)
		{
			sum += fabs(x[i]);
		}
		i = m;
		while(i < size)
		{
			sum += fabs(x[i]) + fabs(x[i + 1]) + fabs(x[i + 2]) + fabs(x[i + 3]) + fabs(x[i + 4]) + fabs(x[i + 5]);
			i += step;
		}
	}
	else
	{
		unsigned int x_index = 0;
		for(unsigned int i = 0 ; i < size ; i++)
		{
			sum += fabs(x[x_index]);
			x_index += x_inc;
		}
	}
	return sum;
}
double ddot(unsigned int size,unsigned int x_inc,const double* x,unsigned int y_inc,const double* y)
{
	// This function computes the dot product of the arrays X and Y such that the dot product is the 
	// sum of the products of the corresponding array elements spaced x_inc and y_inc elements apart and 
	// the sum runs up to N product pairs. In case of x_inc = y_inc = 1, this is the regular dot product 
	// between two arrays. The specification of x_inc and y_inc is useful so that this function can be 
	// used to compute the dot product between a row and a column vector for example in two matrices. 
	// In this case, x_inc and y_inc can be used to specify the column elements if its matrix is stored 
	// row wise. 
	double product = 0.0;
	if(size == 0)					return product;
	if((x_inc == 1) && (y_inc == 1))
	{
		// unroll the loops 5 at a time
		unsigned int step = 5;
		unsigned int i = 0;
		unsigned int m = size%step;
		for(i = 0 ; i < m ; i++)
		{
			product += x[i]*y[i];
		}
		i = m;
		while(i < size)
		{
			product += x[i]*y[i] + x[i + 1]*y[i + 1] + x[i + 2]*y[i + 2] + x[i + 3]*y[i + 3] + x[i + 4]*y[i + 4];
			i += step;
		}
	}
	else
	{
		unsigned int x_index = 0;
		unsigned int y_index = 0;
		for(unsigned int i = 0 ; i < size ; i++)
		{
			product += x[x_index]*y[y_index];
			x_index += x_inc;
			y_index += y_inc;
		}
	}
	return product;
}
double dnrm2(unsigned int size,unsigned int x_inc,const double* x)
{
	// This function computes the square of the Euclidean norm of N values of the array X 
	// entries spaced x_inc apart.
	double norm = 0.0;
	if(size == 0)					return norm;
	if(x_inc == 0)					return norm;
	if(size == 1)					return (x[0]*x[0]);
	unsigned int x_index = 0;
	for(unsigned int i = 0 ; i < size ; i++)
	{
		norm += x[x_index]*x[x_index];
		x_index += x_inc;
	}
	return norm;
}
double dnrm(unsigned int size,unsigned int x_inc,const double* x)
{
	// This function computes the Euclidean norm of N values of the array X 
	// entries spaced x_inc apart.
	return sqrt(dnrm2(size,x_inc,x));
}
unsigned int idamax(unsigned int size,unsigned int x_inc,const double* x)
{
	// This function finds the index of the first element in the array X that has the maximum 
	// absolute value out of N elements for the elements that are spaced x_inc apart.
	unsigned int index = 0;
	if(size < 2)				return index;
	if(x_inc == 0)				return index;
	double max = 0.0;
	if(x_inc == 1)
	{
		max = fabs(x[0]);
		for(unsigned int i = 1 ; i < size ; i++)
		{
			if(fabs(x[i]) > max)
			{
				max = fabs(x[i]);
				index = i;
			}
		}
	}
	else
	{
		unsigned int x_index = 0;
		max = fabs(x[0]);
		for(unsigned int i = 0 ; i < size ; i++)
		{
			if(fabs(x[x_index]) > max)
			{
				max = fabs(x[x_index]);
				index = i;
			}
			x_index += x_inc;
		}
	}
	return index;
}
void dgemv(unsigned int row_count,unsigned int column_count,double alpha,const double* matrix,unsigned int x_inc,const double* x,double beta,unsigned int y_inc,double* y,int operation)
{
	// This function performs the operation 
	// y = alpha*A*x + beta*y if Operation = 1
	// or 
	// y = alpha*A^T*x + beta*y if Operation = 2
	// in both cases, the results are stored in Y. 
	// A is an MxN matrix, x_inc and y_inc are used to set the entry spacing in the X and Y arrays 
	// respectively, alpha and beta are scalars
	if(x_inc == 0)			return;
	if(y_inc == 0)			return;
	if(row_count == 0)		return;
	if(column_count == 0)	return;
	if(operation < 1)		return;
	if(operation > 2)		return;
	// start by forming the product beta*y if needed and store in Y
	if(operation == 2)		dscal(column_count,y_inc,y,beta);
	else					dscal(row_count,y_inc,y,beta);
	// then form the product alpha*A*x or alpha*A^T*x as instructed if needed and add it to Y
	if(fabs(alpha) > BLAS_tolerance)
	{
		if(operation == 2)
		{
			// form alpha*A^T*x and add it to Y
			if(y_inc == 1)
			{
				for(unsigned int i = 0 ; i < column_count ; i++)
				{
					y[i] += alpha*ddot(row_count,column_count,matrix + i,x_inc,x);
				}
			}
			else
			{
				unsigned int y_index = 0;
				for(unsigned int i = 0 ; i < column_count ; i++)
				{
					y[y_index] += alpha*ddot(row_count,column_count,matrix + i,x_inc,x);
					y_index += y_inc;
				}
			}
		}
		else
		{
			// form alpha*A*x and add it to Y
			if(y_inc == 1)
			{
				for(unsigned int i = 0 ; i < row_count ; i++)
				{
					y[i] += alpha*ddot(column_count,1,matrix + i*column_count,x_inc,x);
				}
			}
			else
			{
				unsigned int y_index = 0;
				for(unsigned int i = 0 ; i < row_count ; i++)
				{
					y[y_index] += alpha*ddot(column_count,1,matrix + i*column_count,x_inc,x);
					y_index += y_inc;
				}
			}
		}
	}
}
void dtrmv(unsigned int size,const double* matrix,unsigned int x_inc,const double* x,unsigned int y_inc,double* y,int matrix_type)
{
	// This function performs the operation 
	// y <-- A*x 
	// or 
	// y <-- A^T*x
	// depending on the matrix type where A is an nxn triangular matrix and X is a vector
	// x_inc and y_inc are used to set the entry spacing in the X and Y arrays
	// matrix types:
	// 1: upper triangular, don't transpose, with a non-unit diagonal
	// 2: upper triangular, don't transpose, with a unit diagonal
	// 3: lower triangular, don't transpose, with a non-unit diagonal
	// 4: lower triangular, don't transpose, with a unit diagonal
	// 5: upper triangular, transpose, with a non-unit diagonal
	// 6: upper triangular, transpose, with a unit diagonal
	// 7: lower triangular, transpose, with a non-unit diagonal
	// 8: lower triangular, transpose, with a unit diagonal
	if(matrix_type < 1)		return;
	if(matrix_type > 8)		return;
	if(x_inc == 0)			return;
	if(y_inc == 0)			return;
	if(size == 0)			return;
	double temp = 0.0;
	// The multiplication order is as follows, instead of performing the regular row-wise 
	// matrix vector product, we will go over all the columns of the A matrix, multiply the jth 
	// column by the jth scalar in the X vector and accumulate the product in the result. 
	// While this does not save us any operations, we are able to skip computations if the jth 
	// vector element is zero, in which case the entire scalar-column product is skipped. This 
	// works only for triangular matrices though because otherwise, the accumulation would pollute 
	// the output array (which is the same as the input array) before future iterations get to 
	// use the original values needed for their execution. 
	// use the matrix as is
	unsigned int x_column_index = 0;
	unsigned int y_column_index = 0;
	unsigned int y_row_index = 0;
	if(matrix_type < 5)
	{
		// use the matrix as is
		if(matrix_type < 3)
		{
			// this is an upper triangular matrix
			if(x_inc == 1)
			{
				for(unsigned int j = 0 ; j < size ; j++)
				{
					if(fabs(x[j]) < BLAS_tolerance)		continue;
					// don't we need to clear y[j] before we start accumulation ? 
					// no actually, its old value, after scaling if needed, will be 
					// used as the initial accumulation value since it is in this 
					// iteration that the accumulation will begin by setting the 
					// initial value (or just leaving it as is) and then adding other 
					// values to it in future iterations.
					temp = x[j];
					y_row_index = 0;
					for(unsigned int i = 0 ; i < j ; i++)
					{
						y[y_row_index] += temp*matrix[i*size + j];
						y_row_index += y_inc;
					}
					// why aren't we incrementing here ? because this is the first 
					// time this vector element is modified, if the main diagonal 
					// is one, then leave it as is because it will be the initial 
					// value for future accumulations which will begin next iteration, 
					// if it the main diagonal element is not one, scale it and use 
					// it as the start value for future accumulation.
					if(matrix_type < 2)		y[y_column_index] = temp*matrix[j*size + j];
					else 					y[y_column_index] = temp;
					y_column_index += y_inc;
				}
			}
			else
			{
				for(unsigned int j = 0 ; j < size ; j++)
				{
					if(fabs(x[x_column_index]) < BLAS_tolerance)		continue;
					temp = x[x_column_index];
					y_row_index = 0;
					for(unsigned int i = 0 ; i < j ; i++)
					{
						y[y_row_index] += temp*matrix[i*size + j];
						y_row_index += y_inc;
					}
					if(matrix_type < 2)		y[y_column_index] = temp*matrix[j*size + j];
					else 					y[y_column_index] = temp;
					x_column_index += x_inc;
					y_column_index += y_inc;
				}
			}
		}
		else
		{
			// this is a lower triangular matrix
			if(x_inc == 1)
			{
				y_column_index = size*y_inc;
				for(unsigned int j = (size - 1) ; j >= 0 ; j--)
				{
					y_column_index -= y_inc;
					if(fabs(x[j]) < BLAS_tolerance)
					{
						if(j == 0)		break;
						continue;
					}
					temp = x[j];
					y_row_index = size*y_inc;
					for(unsigned int i = (size - 1) ; i > j ; i--)
					{
						y_row_index -= y_inc;
						y[y_row_index] += temp*matrix[i*size + j];
					}
					if(matrix_type < 4)		y[y_column_index] = temp*matrix[j*size + j];
					else 					y[y_column_index] = temp;
					if(j == 0)		break;
				}
			}
			else
			{
				y_column_index = size*y_inc;
				x_column_index = size*x_inc;
				for(unsigned int j = (size - 1) ; j >= 0 ; j--)
				{
					y_column_index -= y_inc;
					x_column_index -= x_inc;
					if(fabs(x[x_column_index]) < BLAS_tolerance)
					{
						if(j == 0)		break;
						continue;
					}
					temp = x[x_column_index];
					y_row_index = size*y_inc;
					for(unsigned int i = (size - 1) ; i > j ; i--)
					{
						y_row_index -= y_inc;
						y[y_row_index] += temp*matrix[i*size + j];
					}
					if(matrix_type < 4)		y[y_column_index] = temp*matrix[j*size + j];
					else 					y[y_column_index] = temp;
					if(j == 0)		break;
				}
			}
		}
	}
	else
	{
		// use the matrix transpose
		if(matrix_type < 7)
		{
			// this is an upper triangular matrix
			if(x_inc == 1)
			{
				y_column_index = size*y_inc;
				for(unsigned int j = (size - 1) ; j >= 0 ; j--)
				{
					y_column_index -= y_inc;
					if(fabs(x[j]) < BLAS_tolerance)
					{
						if(j == 0)		break;
						continue;
					}
					temp = x[j];
					y_row_index = size*y_inc;
					for(unsigned int i = (size - 1) ; i > j ; i--)
					{
						y_row_index -= y_inc;
						y[y_row_index] += temp*matrix[j*size + i];
					}
					if(matrix_type < 6)		y[y_column_index] = temp*matrix[j*size + j];
					else 					y[y_column_index] = temp;
					if(j == 0)				break;
				}
			}
			else
			{
				y_column_index = size*y_inc;
				x_column_index = size*x_inc;
				for(unsigned int j = (size - 1) ; j >= 0 ; j--)
				{
					y_column_index -= y_inc;
					x_column_index -= x_inc;
					if(fabs(x[x_column_index]) < BLAS_tolerance)
					{
						if(j == 0)		break;
						continue;
					}
					temp = x[x_column_index];
					y_row_index = size*y_inc;
					for(unsigned int i = (size - 1) ; i > j ; i--)
					{
						y_row_index -= y_inc;
						y[y_row_index] += temp*matrix[j*size + i];
					}
					if(matrix_type < 6)		y[y_column_index] = temp*matrix[j*size + j];
					else 					y[y_column_index] = temp;
					if(j == 0)				break;
				}
			}
		}
		else
		{
			// this is a lower triangular matrix
			if(x_inc == 1)
			{
				y_column_index = 0;
				for(unsigned int j = 0 ; j < size ; j++)
				{
					if(fabs(x[j]) < BLAS_tolerance)		continue;
					temp = x[j];
					y_row_index = 0;
					for(unsigned int i = 0 ; i < j ; i++)
					{
						y[y_row_index] += temp*matrix[j*size + i];
						y_row_index += y_inc;
					}
					if(matrix_type < 8)		y[y_column_index] = temp*matrix[j*size + j];
					else 					y[y_column_index] = temp;
					y_column_index += y_inc;
				}
			}
			else
			{
				y_column_index = 0;
				x_column_index = 0;
				for(unsigned int j = 0 ; j < size ; j++)
				{
					if(fabs(x[x_column_index]) < BLAS_tolerance)		continue;
					temp = x[x_column_index];
					y_row_index = 0;
					for(unsigned int i = 0 ; i < j ; i++)
					{
						y[y_row_index] += temp*matrix[j*size + i];
						y_row_index += y_inc;
					}
					if(matrix_type < 8)		y[y_column_index] = temp*matrix[j*size + j];
					else 					y[y_column_index] = temp;
					x_column_index += x_inc;
					y_column_index += y_inc;
				}
			}
		}
	}
}
void dtrsv(unsigned int size,const double* matrix,unsigned int x_inc,double* x,int matrix_type)
{
	// This function solves the system of equations of the form
	// A*x = b
	// or 
	// A^T*x = b
	// depending on the matrix type where A is an nxn triangular matrix and b is a vector. 
	// The X input array holds the values for the b vector on entry and is overwritten with 
	// the solution on exit. x_inc is used to set the entry spacing in the X array.
	// matrix types:
	// 1: upper triangular, don't transpose, with a non-unit diagonal
	// 2: upper triangular, don't transpose, with a unit diagonal
	// 3: lower triangular, don't transpose, with a non-unit diagonal
	// 4: lower triangular, don't transpose, with a unit diagonal
	// 5: upper triangular, transpose, with a non-unit diagonal
	// 6: upper triangular, transpose, with a unit diagonal
	// 7: lower triangular, transpose, with a non-unit diagonal
	// 8: lower triangular, transpose, with a unit diagonal
	if(matrix_type < 1)		return;
	if(matrix_type > 8)		return;
	if(x_inc == 0)			return;
	if(size == 0)			return;
	unsigned int x_row_index = 0;
	unsigned int x_column_index = 0;
	if(matrix_type < 5)
	{
		// use the matrix as is
		if(matrix_type < 3)
		{
			// this is an upper triangular matrix
			if(x_inc == 1)
			{
				for(int j = (size - 1) ; j >= 0 ; j--)
				{
					if(fabs(x[j]) < BLAS_tolerance)		continue;
					if(matrix_type < 2)			x[j] = x[j]/matrix[j*size + j];
					for(int i = (j - 1) ; i >= 0 ; i--)
					{
						x[i] -= matrix[i*size + j]*x[j];
					}
				}
			}
			else
			{
				x_row_index = 0;
				x_column_index = size*x_inc;
				for(int j = (size - 1) ; j >= 0 ; j--)
				{
					x_column_index -= x_inc;
					if(fabs(x[x_column_index]) < BLAS_tolerance)		continue;
					if(matrix_type < 2)			x[x_column_index] = x[x_column_index]/matrix[j*size + j];
					x_row_index = j*x_inc;
					for(int i = (j - 1) ; i >= 0 ; i--)
					{
						x_row_index -= x_inc;
						x[x_row_index] -= matrix[i*size + j]*x[x_column_index];
					}
					
				}
			}
		}
		else
		{
			// this is a lower triangular matrix
			if(x_inc == 1)
			{
				for(unsigned int j = 0 ; j < size ; j++)
				{
					if(fabs(x[j]) < BLAS_tolerance)		continue;
					if(matrix_type < 4)			x[j] = x[j]/matrix[j*size + j];
					for(unsigned int i = (j + 1) ; i < size ; i++)
					{
						x[i] -= matrix[i*size + j]*x[j];
					}
				}
			}
			else
			{
				x_row_index = 0;
				x_column_index = 0;
				for(unsigned int j = 0 ; j < size ; j++)
				{
					if(fabs(x[x_column_index]) < BLAS_tolerance)	continue;
					if(matrix_type < 4)			x[x_column_index] = x[x_column_index]/matrix[j*size + j];
					x_row_index = j*x_inc;
					for(unsigned int i = (j + 1) ; i < size ; i++)
					{
						x_row_index += x_inc;
						x[x_row_index] -= matrix[i*size + j]*x[x_column_index];
					}
					x_column_index += x_inc;
				}
			}
		}
	}
	else
	{
		// use the matrix transpose
		if(matrix_type < 7)
		{
			// this is an upper triangular matrix
			if(x_inc == 1)
			{
				for(unsigned int j = 0 ; j < size ; j++)
				{
					if(fabs(x[j]) < BLAS_tolerance)			continue;
					if(matrix_type < 6)			x[j] = x[j]/matrix[j*size + j];
					for(unsigned int i = (j + 1) ; i < size ; i++)
					{
						x[i] -= matrix[j*size + i]*x[j];
					}
				}
			}
			else
			{
				x_row_index = 0;
				x_column_index = 0;
				for(unsigned int j = 0 ; j < size ; j++)
				{
					if(fabs(x[x_column_index]) < BLAS_tolerance)		continue;
					if(matrix_type < 6)			x[x_column_index] = x[x_column_index]/matrix[j*size + j];
					x_row_index = j*x_inc;
					for(unsigned int i = (j + 1) ; i < size ; i++)
					{
						x_row_index += x_inc;
						x[x_row_index] -= matrix[j*size + i]*x[x_column_index];
					}
					x_column_index += x_inc;
				}
			}
		}
		else
		{
			// this is a lower triangular matrix
			if(x_inc == 1)
			{
				for(int j = (size - 1) ; j >= 0 ; j--)
				{
					if(fabs(x[j]) < BLAS_tolerance)		continue;
					if(matrix_type < 8)			x[j] = x[j]/matrix[j*size + j];
					for(int i = (j - 1) ; i >= 0 ; i--)
					{
						x[i] -= matrix[j*size + i]*x[j];
					}
				}
			}
			else
			{
				x_row_index = 0;
				x_column_index = size*x_inc;
				for(int j = (size - 1) ; j >= 0 ; j--)
				{
					x_column_index -= x_inc;
					if(fabs(x[x_column_index]) < BLAS_tolerance)		continue;
					if(matrix_type < 8)			x[x_column_index] = x[x_column_index]/matrix[j*size + j];
					x_row_index = j*x_inc;
					for(int i = (j - 1) ; i >= 0 ; i--)
					{
						x_row_index -= x_inc;
						x[x_row_index] -= matrix[j*size + i]*x[x_column_index];
					}
				}
			}
		}
	}
}
void ddgmv(unsigned int size,const double* matrix,unsigned int x_inc,const double* x,unsigned int y_inc,double* y,int matrix_type)
{
	// This function multiplies the diagonal matrix A with the column vector X and stores the 
	// result in the column vector Y where A is an nxn triangular matrix and X is a vector
	// x_inc and y_inc are used to set the entry spacing in the X and Y arrays. If the matrix 
	// type is 1, then the full matrix is passed in A, if it is 2, then A is just an array of 
	// the diagonal entries of the matrix A
	if(x_inc == 0)			return;
	if(y_inc == 0)			return;
	if(size == 0)			return;
	if(matrix_type < 1)		return;
	if(matrix_type > 2)		return;
	unsigned int y_column_index = 0;
	unsigned int x_column_index = 0;
	if(matrix_type == 1)
	{
		// A is the full matrix
		if(x_inc == 1)
		{
			for(unsigned int j = 0 ; j < size ; j++)
			{
				y[j] = x[j]*matrix[j*size + j];
			}
		}
		else
		{
			y_column_index = 0;
			x_column_index = 0;
			for(unsigned int j = 0 ; j < size ; j++)
			{
				y[y_column_index] = x[x_column_index]*matrix[j*size + j];
				x_column_index += x_inc;
				y_column_index += y_inc;
			}
		}
	}
	else if(matrix_type == 2)
	{
		// A is the array of diagonal entries
		if(x_inc == 1)
		{
			for(unsigned int j = 0 ; j < size ; j++)
			{
				y[j] = x[j]*matrix[j];
			}
		}
		else
		{
			y_column_index = 0;
			x_column_index = 0;
			for(unsigned int j = 0 ; j < size ; j++)
			{
				y[y_column_index] = x[x_column_index]*matrix[j];
				x_column_index += x_inc;
				y_column_index += y_inc;
			}
		}
	}
}
int test_BLAS(unsigned int size,int comprehensive,double tolerance,unsigned int paranoia)
{
	int pass = 1;
	if(comprehensive)
	{
		if(!test_dcopy(size,paranoia,tolerance))
		{
			printf("BLAS dcopy failed test\n");
			pass = 0;
		}
		if(!test_dswap(size,paranoia,tolerance))
		{
			printf("BLAS dswap failed test\n");
			pass = 0;
		}
		if(!test_dscal(size,paranoia,tolerance))
		{
			printf("BLAS dscal failed test\n");
			pass = 0;
		}
		if(!test_daxpy(size,paranoia,tolerance))
		{
			printf("BLAS daxpy failed test\n");
			pass = 0;
		}
		if(!test_drotg(paranoia,tolerance))
		{
			printf("BLAS drotg failed test\n");
			pass = 0;
		}
		if(!test_drot(paranoia,tolerance))
		{
			printf("BLAS drot failed test\n");
			pass = 0;
		}
		if(!test_dasum())
		{
			printf("BLAS dasum failed test\n");
			pass = 0;
		}
		if(!test_ddot(size,paranoia,tolerance))
		{
			printf("BLAS ddot failed test\n");
			pass = 0;
		}
		if(!test_dnrm2(size,paranoia,tolerance))
		{
			printf("BLAS drnm2 failed test\n");
			pass = 0;
		}
		if(!test_dnrm(size,paranoia,tolerance))
		{
			printf("BLAS dnrm failed test\n");
			pass = 0;
		}
		if(!test_idamax())
		{
			printf("BLAS idamax failed test\n");
			pass = 0;
		}
		if(!test_dgemv(size,paranoia,tolerance))
		{
			printf("BLAS dgemv failed test\n");
			pass = 0;
		}
		if(!test_dtrmv(size,paranoia,tolerance))
		{
			printf("BLAS dtrmv failed test\n");
			pass = 0;
		}
		if(!test_dtrsv(size,paranoia,tolerance))
		{
			printf("BLAS dtrsv failed test\n");
			pass = 0;
		}
		if(!test_ddgmv(size,paranoia,tolerance))
		{
			printf("BLAS ddgmv failed test\n");
			pass = 0;
		}
	}
	else
	{
		if(!test_ddgmv(size,paranoia,tolerance))
		{
			printf("BLAS ddgmv failed test\n");
			pass = 0;
		}
	}
	return pass;
}
int test_dcopy(unsigned int test_size,unsigned int test_paranoia,double tolerance)
{
	int pass = 1;
	// test unit increment copy in source and destination
	unsigned int source_stride = 1;
	unsigned int destination_stride = 1;
	double* source = (double*)malloc(test_size*source_stride*sizeof(double));
	double* destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*source_stride,source);
		dcopy(test_size,source_stride,source,destination_stride,destination);
		pass = compare_vectors(test_size,source_stride,source,destination_stride,destination,1.0,0,0,tolerance);
		if(!pass)			break;
	}
	free(source);
	free(destination);
	if(!pass)				return 0;
	// test unit increment copy in source but not in destination, generate a random number 
	// for strides but keep it under 10.
	source_stride = 1;
	while(destination_stride <= 1)
	{
		destination_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
	}
	source = (double*)malloc(test_size*source_stride*sizeof(double));
	destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*source_stride,source);
		dcopy(test_size,source_stride,source,destination_stride,destination);
		pass = compare_vectors(test_size,source_stride,source,destination_stride,destination,1.0,0,0,tolerance);
		if(!pass)			break;
	}
	free(source);
	free(destination);
	if(!pass)				return 0;
	// test strided increment copy in source but unit increment in destination, generate a random number 
	// for strides but keep it under 10.
	while(source_stride <= 1)
	{
		source_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	destination_stride = 1;
	source = (double*)malloc(test_size*source_stride*sizeof(double));
	destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*source_stride,source);
		dcopy(test_size,source_stride,source,destination_stride,destination);
		pass = compare_vectors(test_size,source_stride,source,destination_stride,destination,1.0,0,0,tolerance);
		if(!pass)			break;
	}
	free(source);
	free(destination);
	if(!pass)				return 0;
	// test strided increment copy in both source and destination, generate a random number 
	// for strides but keep it under 10.
	source_stride = 1;
	destination_stride = 1;
	while(source_stride <= 1)
	{
		source_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	while(destination_stride <= 1)
	{
		destination_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
	}
	source = (double*)malloc(test_size*source_stride*sizeof(double));
	destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*source_stride,source);
		dcopy(test_size,source_stride,source,destination_stride,destination);
		pass = compare_vectors(test_size,source_stride,source,destination_stride,destination,1.0,0,0,tolerance);
		if(!pass)			break;
	}
	free(source);
	free(destination);
	return pass;
}
int test_dswap(unsigned int test_size,unsigned int test_paranoia,double tolerance)
{
	// This function uses dcopy(). It is assumed that it passed the test. 
	int pass = 1;
	// test unit increment swap in source and destination
	unsigned int source_stride = 1;
	unsigned int destination_stride = 1;
	double* source = (double*)malloc(test_size*source_stride*sizeof(double));
	double* original_source = (double*)malloc(test_size*source_stride*sizeof(double));
	double* destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	double* original_destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*source_stride,source);
		dcopy(test_size,source_stride,source,source_stride,original_source);
		vector_random_populate(test_size*destination_stride,destination);
		dcopy(test_size,destination_stride,destination,destination_stride,original_destination);
		dswap(test_size,source_stride,source,destination_stride,destination);
		pass = compare_vectors(test_size,source_stride,source,destination_stride,original_destination,1.0,0,0,tolerance);
		pass = compare_vectors(test_size,destination_stride,destination,source_stride,original_source,1.0,0,0,tolerance);
		if(!pass)			break;
	}
	free(source);
	free(original_source);
	free(destination);
	free(original_destination);
	if(!pass)				return 0;
	// test unit increment swap in source but not in destination, generate a random number 
	// for strides but keep it under 10.
	source_stride = 1;
	while(destination_stride <= 1)
	{
		destination_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
	}
	source = (double*)malloc(test_size*source_stride*sizeof(double));
	original_source = (double*)malloc(test_size*source_stride*sizeof(double));
	destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	original_destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*source_stride,source);
		dcopy(test_size,source_stride,source,source_stride,original_source);
		vector_random_populate(test_size*destination_stride,destination);
		dcopy(test_size,destination_stride,destination,destination_stride,original_destination);
		dswap(test_size,source_stride,source,destination_stride,destination);
		pass = compare_vectors(test_size,source_stride,source,destination_stride,original_destination,1.0,0,0,tolerance);
		pass = compare_vectors(test_size,destination_stride,destination,source_stride,original_source,1.0,0,0,tolerance);
		if(!pass)			break;
	}
	free(source);
	free(original_source);
	free(destination);
	free(original_destination);
	if(!pass)				return 0;
	// test strided increment swap in source but unit increment in destination, generate a random number 
	// for strides but keep it under 10.
	while(source_stride <= 1)
	{
		source_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	destination_stride = 1;
	source = (double*)malloc(test_size*source_stride*sizeof(double));
	original_source = (double*)malloc(test_size*source_stride*sizeof(double));
	destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	original_destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*source_stride,source);
		dcopy(test_size,source_stride,source,source_stride,original_source);
		vector_random_populate(test_size*destination_stride,destination);
		dcopy(test_size,destination_stride,destination,destination_stride,original_destination);
		dswap(test_size,source_stride,source,destination_stride,destination);
		pass = compare_vectors(test_size,source_stride,source,destination_stride,original_destination,1.0,0,0,tolerance);
		pass = compare_vectors(test_size,destination_stride,destination,source_stride,original_source,1.0,0,0,tolerance);
		if(!pass)			break;
	}
	free(source);
	free(original_source);
	free(destination);
	free(original_destination);
	if(!pass)				return 0;
	// test strided increment swap in both source and destination, generate a random number 
	// for strides but keep it under 10.
	source_stride = 1;
	destination_stride = 1;
	while(source_stride <= 1)
	{
		source_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	while(destination_stride <= 1)
	{
		destination_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
	}
	source = (double*)malloc(test_size*source_stride*sizeof(double));
	original_source = (double*)malloc(test_size*source_stride*sizeof(double));
	destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	original_destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*source_stride,source);
		dcopy(test_size,source_stride,source,source_stride,original_source);
		vector_random_populate(test_size*destination_stride,destination);
		dcopy(test_size,destination_stride,destination,destination_stride,original_destination);
		dswap(test_size,source_stride,source,destination_stride,destination);
		pass = compare_vectors(test_size,source_stride,source,destination_stride,original_destination,1.0,0,0,tolerance);
		pass = compare_vectors(test_size,destination_stride,destination,source_stride,original_source,1.0,0,0,tolerance);
		if(!pass)			break;
	}
	free(source);
	free(original_source);
	free(destination);
	free(original_destination);
	return pass;
}
int test_dscal(unsigned int test_size,unsigned int test_paranoia,double tolerance)
{
	// This function uses dcopy(). It is assumed that it passed the test. 
	double factor = 1.0;
	int pass = 1;
	// factors will range from -100.0 to 100.0
	// test unit increment scaling
	unsigned int stride = 1;
	double* working_array = (double*)malloc(test_size*stride*sizeof(double));
	double* original_array = (double*)malloc(test_size*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*stride,working_array);
		factor = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		dcopy(test_size,stride,working_array,1,original_array);
		dscal(test_size,stride,working_array,factor);
		pass = compare_vectors(test_size,stride,working_array,1,original_array,factor,0,0,tolerance);
		if(!pass)
		{
			printf("BLAS Scale test failed for factor %e\n",factor);
			break;
		}
	}
	free(working_array);
	if(!pass)
	{
		free(original_array);
		return 0;
	}
	// test strided increment scaling, generate a random number for strides but keep it
	// under 10.
	while(stride <= 1)
	{
		stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	working_array = (double*)malloc(test_size*stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*stride,working_array);
		factor = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		dcopy(test_size,stride,working_array,1,original_array);
		dscal(test_size,stride,working_array,factor);
		pass = compare_vectors(test_size,stride,working_array,1,original_array,factor,0,0,tolerance);
		if(!pass)
		{
			printf("dscal test failed for factor %e\n",factor);
			break;
		}
	}
	free(working_array);
	free(original_array);
	return pass;
}
int test_daxpy(unsigned int test_size,unsigned int test_paranoia,double tolerance)
{
	// This function uses dcopy(). It is assumed that it passed the test. 
	double factor = 0.0;
	int pass = 1;
	// factors will range from -100.0 to 100.0
	// test unit increment scale and add in source and destination
	unsigned int source_stride = 1;
	unsigned int destination_stride = 1;
	double* source = (double*)malloc(test_size*source_stride*sizeof(double));
	double* original_source = (double*)malloc(test_size*source_stride*sizeof(double));
	double* destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	double* original_destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*source_stride,source);
		dcopy(test_size,source_stride,source,source_stride,original_source);
		vector_random_populate(test_size*destination_stride,destination);
		dcopy(test_size,destination_stride,destination,destination_stride,original_destination);
		factor = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		daxpy(test_size,source_stride,source,factor,destination_stride,destination);
		pass = compare_vectors(test_size,destination_stride,destination,source_stride,original_source,factor,destination_stride,original_destination,tolerance);
		if(!pass)			break;
	}
	free(source);
	free(original_source);
	free(destination);
	free(original_destination);
	if(!pass)				return 0;
	// test unit increment scale and add in source but not in destination, generate a random number 
	// for strides but keep it under 10.
	source_stride = 1;
	while(destination_stride <= 1)
	{
		destination_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
	}
	source = (double*)malloc(test_size*source_stride*sizeof(double));
	original_source = (double*)malloc(test_size*source_stride*sizeof(double));
	destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	original_destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*source_stride,source);
		dcopy(test_size,source_stride,source,source_stride,original_source);
		vector_random_populate(test_size*destination_stride,destination);
		dcopy(test_size,destination_stride,destination,destination_stride,original_destination);
		factor = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		daxpy(test_size,source_stride,source,factor,destination_stride,destination);
		pass = compare_vectors(test_size,destination_stride,destination,source_stride,original_source,factor,destination_stride,original_destination,tolerance);
		if(!pass)			break;
	}
	free(source);
	free(original_source);
	free(destination);
	free(original_destination);
	if(!pass)				return 0;
	// test strided increment scale and add in source but unit increment in destination, generate a random number 
	// for strides but keep it under 10.
	while(source_stride <= 1)
	{
		source_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	destination_stride = 1;
	source = (double*)malloc(test_size*source_stride*sizeof(double));
	original_source = (double*)malloc(test_size*source_stride*sizeof(double));
	destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	original_destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*source_stride,source);
		dcopy(test_size,source_stride,source,source_stride,original_source);
		vector_random_populate(test_size*destination_stride,destination);
		dcopy(test_size,destination_stride,destination,destination_stride,original_destination);
		factor = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		daxpy(test_size,source_stride,source,factor,destination_stride,destination);
		pass = compare_vectors(test_size,destination_stride,destination,source_stride,original_source,factor,destination_stride,original_destination,tolerance);
		if(!pass)			break;
	}
	free(source);
	free(original_source);
	free(destination);
	free(original_destination);
	if(!pass)				return 0;
	// test strided increment scale and add in both source and destination, generate a random number 
	// for strides but keep it under 10.
	source_stride = 1;
	destination_stride = 1;
	while(source_stride <= 1)
	{
		source_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	while(destination_stride <= 1)
	{
		destination_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
	}
	source = (double*)malloc(test_size*source_stride*sizeof(double));
	original_source = (double*)malloc(test_size*source_stride*sizeof(double));
	destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	original_destination = (double*)malloc(test_size*destination_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*source_stride,source);
		dcopy(test_size,source_stride,source,source_stride,original_source);
		vector_random_populate(test_size*destination_stride,destination);
		dcopy(test_size,destination_stride,destination,destination_stride,original_destination);
		factor = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		daxpy(test_size,source_stride,source,factor,destination_stride,destination);
		pass = compare_vectors(test_size,destination_stride,destination,source_stride,original_source,factor,destination_stride,original_destination,tolerance);
		if(!pass)			break;
	}
	free(source);
	free(original_source);
	free(destination);
	free(original_destination);
	return pass;
}
int test_drotg(unsigned int test_paranoia,double tolerance)
{
	int pass = 1;
	double a = 0.0;
	double b = 0.0;
	double angle_cosine = 0.0;
	double angle_sine = 0.0;
	double x = 0.0;
	double y = 0.0;
	double r = 0.0;
	// generate A and B in the range -100.0 to 100.0
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		a = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		b = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		r = sqrt(a*a + b*b);
		drotg(a,b,&angle_cosine,&angle_sine);
		x = angle_cosine*a + angle_sine*b;
		y = -angle_sine*a + angle_cosine*b;
		if(fabs(y) > tolerance)					pass = 0;
		if(fabs(fabs(x) - r) > tolerance)		pass = 0;
		r = angle_cosine*angle_cosine + angle_sine*angle_sine;
		if(fabs(r - 1.0) > tolerance)			pass = 0;
		if(!pass)									break;
	}
	return pass;
}
int test_drot(unsigned int test_paranoia,double tolerance)
{
	// drot() passes if drotg() passes
	return test_drotg(test_paranoia,tolerance);
}
int test_dasum()
{
	// The test for dasum will be another implementation of dasum
	return 1;
}
int test_ddot(unsigned int test_size,unsigned int test_paranoia,double tolerance)
{
	double product = 0.0;
	int pass = 1;
	// test unit increment dot product in X and Y
	unsigned int x_stride = 1;
	unsigned int y_stride = 1;
	double* x = (double*)malloc(test_size*x_stride*sizeof(double));
	double* y = (double*)malloc(test_size*y_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*x_stride,x);
		vector_random_populate(test_size*y_stride,y);
		product = ddot(test_size,x_stride,x,y_stride,y);
		for(unsigned int j = 0 ; j < test_size ; j++)
		{
			product -= x[j*x_stride]*y[j*y_stride];
		}
		pass = (fabs(product) < tolerance);
		if(!pass)
		{
			printf("BLAS ddot test failed with error : %e\n",product);
			break;
		}
	}
	free(x);
	free(y);
	if(!pass)				return 0;
	// test unit increment dot product in X but not in Y, generate a random number 
	// for strides but keep it under 10.
	x_stride = 1;
	while(y_stride <= 1)
	{
		y_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
	}
	x = (double*)malloc(test_size*x_stride*sizeof(double));
	y = (double*)malloc(test_size*y_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*x_stride,x);
		vector_random_populate(test_size*y_stride,y);
		product = ddot(test_size,x_stride,x,y_stride,y);
		for(unsigned int j = 0 ; j < test_size ; j++)
		{
			product -= x[j*x_stride]*y[j*y_stride];
		}
		pass = (fabs(product) < tolerance);
		if(!pass)
		{
			printf("BLAS ddot test failed with error : %e\n",product);
			break;
		}
	}
	free(x);
	free(y);
	if(!pass)				return 0;
	// test strided increment dot product in X but unit increment in Y, generate a random number 
	// for strides but keep it under 10.
	while(x_stride <= 1)
	{
		x_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	y_stride = 1;
	x = (double*)malloc(test_size*x_stride*sizeof(double));
	y = (double*)malloc(test_size*y_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*x_stride,x);
		vector_random_populate(test_size*y_stride,y);
		product = ddot(test_size,x_stride,x,y_stride,y);
		for(unsigned int j = 0 ; j < test_size ; j++)
		{
			product -= x[j*x_stride]*y[j*y_stride];
		}
		pass = (fabs(product) < tolerance);
		if(!pass)
		{
			printf("BLAS ddot test failed with error : %e\n",product);
			break;
		}
	}
	free(x);
	free(y);
	if(!pass)				return 0;
	// test strided increment dot product in both X and Y, generate a random number 
	// for strides but keep it under 10.
	x_stride = 1;
	y_stride = 1;
	while(x_stride <= 1)
	{
		x_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	while(y_stride <= 1)
	{
		y_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
	}
	x = (double*)malloc(test_size*x_stride*sizeof(double));
	y = (double*)malloc(test_size*y_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*x_stride,x);
		vector_random_populate(test_size*y_stride,y);
		product = ddot(test_size,x_stride,x,y_stride,y);
		for(unsigned int j = 0 ; j < test_size ; j++)
		{
			product -= x[j*x_stride]*y[j*y_stride];
		}
		pass = (fabs(product) < tolerance);
		if(!pass)
		{
			printf("BLAS ddot test failed with error : %e\n",product);
			break;
		}
	}
	free(x);
	free(y);
	return pass;
}
int test_dnrm2(unsigned int test_size,unsigned int test_paranoia,double tolerance)
{
	// This function uses ddot(). It is assumed that it passed the test. 
	double squared_norm = 0.0;
	double product = 0.0;
	int pass = 1;
	// test unit increment array squared norm
	unsigned int stride = 1;
	double* x = (double*)malloc(test_size*stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*stride,x);
		squared_norm = dnrm2(test_size,stride,x);
		product = ddot(test_size,stride,x,stride,x);
		pass = fabs(squared_norm - product) < tolerance;
		if(!pass)				break;
	}
	free(x);
	if(!pass)				return 0;
	// test strided array squared norm
	stride = 1;
	while(stride <= 1)
	{
		stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
	}
	x = (double*)malloc(test_size*stride*sizeof(double)); 
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*stride,x);
		squared_norm = dnrm2(test_size,stride,x);
		product = ddot(test_size,stride,x,stride,x);
		pass = fabs(squared_norm - product) < tolerance;
		if(!pass)				break;
	}
	free(x);
	return pass;
}
int test_dnrm(unsigned int test_size,unsigned int test_paranoia,double tolerance)
{
	// This function uses dnrm2(). It is assumed that it passed the test. 
	double squared_norm = 0.0;
	double norm = 0.0;
	int pass = 1;
	// test unit increment array squared norm
	unsigned int stride = 1;
	double* x = (double*)malloc(test_size*stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*stride,x);
		squared_norm = dnrm2(test_size,stride,x);
		norm = dnrm(test_size,stride,x);
		pass = fabs(squared_norm - norm*norm) < tolerance;
		if(!pass)				break;
	}
	free(x);
	if(!pass)				return 0;
	// test strided array squared norm
	stride = 1;
	while(stride <= 1)
	{
		stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
	}
	x = (double*)malloc(test_size*stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		vector_random_populate(test_size*stride,x);
		squared_norm = dnrm2(test_size,stride,x);
		norm = dnrm(test_size,stride,x);
		pass = fabs(squared_norm - norm*norm) < tolerance;
		if(!pass)				break;
	}
	free(x);
	return pass;
}
int test_idamax()
{
	// The test for MaxAbsoluteIndex will be another implementation of MaxAbsoluteIndex
	return 1;
}
int test_dgemv(unsigned int test_size,unsigned int test_paranoia,double tolerance)
{
	// This function uses dcopy() and ddot(). It is assumed that
	// they passed the test. 
	double alpha = 0.0;
	double beta = 0.0;
	double target = 0.0;
	int pass = 1;
	// alpha and beta factors will be chosen in the range -100.0 to 100.0
	// test unit increment matrix vector product in X and Y
	unsigned int x_stride = 1;
	unsigned int y_stride = 1;
	double* matrix = (double*)malloc(test_size*test_size*sizeof(double));
	double* x = (double*)malloc(test_size*x_stride*sizeof(double));
	double* y = (double*)malloc(test_size*y_stride*sizeof(double));
	double* original_y = (double*)malloc(test_size*y_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		matrix_random_populate(test_size,test_size,matrix,0);
		vector_random_populate(test_size*x_stride,x);
		vector_random_populate(test_size*y_stride,y);
		dcopy(test_size,y_stride,y,y_stride,original_y);
		alpha = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		beta = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		// run both no-transpose and transpose tests
		dgemv(test_size,test_size,alpha,matrix,x_stride,x,beta,y_stride,y,1);
		for(unsigned int j = 0 ; j < test_size ; j++)
		{
			target = alpha*ddot(test_size,1,matrix + j*test_size,x_stride,x) + beta*original_y[j*y_stride];
			if(fabs(target - y[j*y_stride]) > tolerance)
			{
				pass = 0;
				break;
			}
		}
		if(!pass)			break;
		// generate a new Y vector because the old one was overwritten
		vector_random_populate(test_size*y_stride,y);
		dcopy(test_size,y_stride,y,y_stride,original_y);
		dgemv(test_size,test_size,alpha,matrix,x_stride,x,beta,y_stride,y,2);
		for(unsigned int j = 0 ; j < test_size ; j++)
		{
			target = alpha*ddot(test_size,test_size,matrix + j,x_stride,x) + beta*original_y[j*y_stride];
			if(fabs(target - y[j*y_stride]) > tolerance)
			{
				pass = 0;
				break;
			}
		}
		if(!pass)			break;
	}
	free(x);
	free(y);
	free(original_y);
	if(!pass)
	{
		free(matrix);
		return 0;
	}
	// test unit increment matrix vector product in X but not in Y, generate a random number 
	// for strides but keep it under 10.
	x_stride = 1;
	while(y_stride <= 1)
	{
		y_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
	}
	x = (double*)malloc(test_size*x_stride*sizeof(double));
	y = (double*)malloc(test_size*y_stride*sizeof(double));
	original_y = (double*)malloc(test_size*y_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		matrix_random_populate(test_size,test_size,matrix,0);
		vector_random_populate(test_size*x_stride,x);
		vector_random_populate(test_size*y_stride,y);
		dcopy(test_size,y_stride,y,y_stride,original_y);
		alpha = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		beta = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		// run both no-transpose and transpose tests
		dgemv(test_size,test_size,alpha,matrix,x_stride,x,beta,y_stride,y,1);
		for(unsigned int j = 0 ; j < test_size ; j++)
		{
			target = alpha*ddot(test_size,1,matrix + j*test_size,x_stride,x) + beta*original_y[j*y_stride];
			if(fabs(target - y[j*y_stride]) > tolerance)
			{
				pass = 0;
				break;
			}
		}
		if(!pass)			break;
		// generate a new Y vector because the old one was overwritten
		vector_random_populate(test_size*y_stride,y);
		dcopy(test_size,y_stride,y,y_stride,original_y);
		dgemv(test_size,test_size,alpha,matrix,x_stride,x,beta,y_stride,y,2);
		for(unsigned int j = 0 ; j < test_size ; j++)
		{
			target = alpha*ddot(test_size,test_size,matrix + j,x_stride,x) + beta*original_y[j*y_stride];
			if(fabs(target - y[j*y_stride]) > tolerance)
			{
				pass = 0;
				break;
			}
		}
		if(!pass)			break;
	}
	free(x);
	free(y);
	free(original_y);
	if(!pass)
	{
		free(matrix);
		return 0;
	}
	// test strided increment matrix vector product in X but unit increment in Y, generate a random number 
	// for strides but keep it under 10.
	while(x_stride <= 1)
	{
		x_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	y_stride = 1;
	x = (double*)malloc(test_size*x_stride*sizeof(double));
	y = (double*)malloc(test_size*y_stride*sizeof(double));
	original_y = (double*)malloc(test_size*y_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		matrix_random_populate(test_size,test_size,matrix,0);
		vector_random_populate(test_size*x_stride,x);
		vector_random_populate(test_size*y_stride,y);
		dcopy(test_size,y_stride,y,y_stride,original_y);
		alpha = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		beta = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		// run both no-transpose and transpose tests
		dgemv(test_size,test_size,alpha,matrix,x_stride,x,beta,y_stride,y,1);
		for(unsigned int j = 0 ; j < test_size ; j++)
		{
			target = alpha*ddot(test_size,1,matrix + j*test_size,x_stride,x) + beta*original_y[j*y_stride];
			if(fabs(target - y[j*y_stride]) > tolerance)
			{
				pass = 0;
				break;
			}
		}
		if(!pass)			break;
		// generate a new Y vector because the old one was overwritten
		vector_random_populate(test_size*y_stride,y);
		dcopy(test_size,y_stride,y,y_stride,original_y);
		dgemv(test_size,test_size,alpha,matrix,x_stride,x,beta,y_stride,y,2);
		for(unsigned int j = 0 ; j < test_size ; j++)
		{
			target = alpha*ddot(test_size,test_size,matrix + j,x_stride,x) + beta*original_y[j*y_stride];
			if(fabs(target - y[j*y_stride]) > tolerance)
			{
				pass = 0;
				break;
			}
		}
		if(!pass)			break;
	}
	free(x);
	free(y);
	free(original_y);
	if(!pass)
	{
		free(matrix);
		return 0;
	}
	// test strided increment copy in both source and destination, generate a random number 
	// for strides but keep it under 10.
	x_stride = 1;
	y_stride = 1;
	while(x_stride <= 1)
	{
		x_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	while(y_stride <= 1)
	{
		y_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
	}
	x = (double*)malloc(test_size*x_stride*sizeof(double));
	y = (double*)malloc(test_size*y_stride*sizeof(double));
	original_y = (double*)malloc(test_size*y_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		matrix_random_populate(test_size,test_size,matrix,0);
		vector_random_populate(test_size*x_stride,x);
		vector_random_populate(test_size*y_stride,y);
		dcopy(test_size,y_stride,y,y_stride,original_y);
		alpha = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		beta = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
		// run both no-transpose and transpose tests
		dgemv(test_size,test_size,alpha,matrix,x_stride,x,beta,y_stride,y,1);
		for(unsigned int j = 0 ; j < test_size ; j++)
		{
			target = alpha*ddot(test_size,1,matrix + j*test_size,x_stride,x) + beta*original_y[j*y_stride];
			if(fabs(target - y[j*y_stride]) > tolerance)
			{
				pass = 0;
				break;
			}
		}
		if(!pass)			break;
		// generate a new Y vector because the old one was overwritten
		vector_random_populate(test_size*y_stride,y);
		dcopy(test_size,y_stride,y,y_stride,original_y);
		dgemv(test_size,test_size,alpha,matrix,x_stride,x,beta,y_stride,y,2);
		for(unsigned int j = 0 ; j < test_size ; j++)
		{
			target = alpha*ddot(test_size,test_size,matrix + j,x_stride,x) + beta*original_y[j*y_stride];
			if(fabs(target - y[j*y_stride]) > tolerance)
			{
				pass = 0;
				break;
			}
		}
		if(!pass)			break;
	}
	free(x);
	free(y);
	free(original_y);
	free(matrix);
	return pass;
}
int test_dtrmv(unsigned int test_size,unsigned int test_paranoia,double tolerance)
{
	// This function uses dgemv(). It is assumed that it passed 
	// the test.
	int all_pass = 1;
	int pass = 1;
	// run tests on all triangular matrix types
	// test unit increment matrix vector product in X
	unsigned int x_stride = 1;
	unsigned int y_stride = 1;
	double* matrix = (double*)malloc(test_size*test_size*sizeof(double));
	double* x = (double*)malloc(test_size*x_stride*sizeof(double));
	double* y = (double*)malloc(test_size*y_stride*sizeof(double));
	double* full_x = (double*)malloc(test_size*x_stride*sizeof(double));
	int matrix_types[8] = {1,5,2,6,3,7,4,8};
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		// go over all matrix types
		for(unsigned int j = 0 ; j < 8 ; j++)
		{
			// XX calculation
			if(matrix_types[j] == 1)
			{
				matrix_random_populate(test_size,test_size,matrix,0);
				make_upper_triangular(test_size,matrix);
			}
			else if(matrix_types[j] == 3)
			{
				matrix_random_populate(test_size,test_size,matrix,0);
				make_lower_triangular(test_size,matrix);
			}
			if((matrix_types[j]%2 == 0))				make_unit_diagonal(test_size,matrix);
			vector_random_populate(test_size*x_stride,x);
			vector_random_populate(test_size*y_stride,y);
			// better set this array to zero to save on scaling during full matrix vector product
			set_value(test_size*x_stride,full_x,0.0);
			dgemv(test_size,test_size,1.0,matrix,x_stride,x,1.0,x_stride,full_x,(matrix_types[j] - 1)/4 + 1);
			dtrmv(test_size,matrix,x_stride,x,y_stride,y,matrix_types[j]);
			pass = compare_vectors(test_size,x_stride,full_x,y_stride,y,1.0,0,0,tolerance);
			if(!pass)			printf("unit stride type %d dtrmv() test failed\n",matrix_types[j]);
			all_pass &= pass;
		}
	}
	free(x);
	free(y);
	free(full_x);
	// test strided increment matrix vector product in X and Y, generate a random number 
	// for strides but keep it under 10.
	while(x_stride <= 1)
	{
		x_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	while(y_stride <= 1)
	{
		y_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	x = (double*)malloc(test_size*x_stride*sizeof(double));
	y = (double*)malloc(test_size*y_stride*sizeof(double));
	full_x = (double*)malloc(test_size*x_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		// go over all matrix types
		for(unsigned int j = 0 ; j < 8 ; j++)
		{
			// XX calculation
			if(matrix_types[j] == 1)
			{
				matrix_random_populate(test_size,test_size,matrix,0);
				make_upper_triangular(test_size,matrix);
			}
			else if(matrix_types[j] == 3)
			{
				matrix_random_populate(test_size,test_size,matrix,0);
				make_lower_triangular(test_size,matrix);
			}
			if((matrix_types[j]%2 == 0))				make_unit_diagonal(test_size,matrix);
			vector_random_populate(test_size*x_stride,x);
			vector_random_populate(test_size*y_stride,y);
			// better set this array to zero to save on scaling during full matrix vector product
			set_value(test_size*x_stride,full_x,0.0);
			dgemv(test_size,test_size,1.0,matrix,x_stride,x,1.0,x_stride,full_x,(matrix_types[j] - 1)/4 + 1);
			dtrmv(test_size,matrix,x_stride,x,y_stride,y,matrix_types[j]);
			pass = compare_vectors(test_size,x_stride,full_x,y_stride,y,1.0,0,0,tolerance);
			if(!pass)			printf("strided type %d dtrmv() test failed\n",matrix_types[j]);
			all_pass &= pass;
		}
	}
	free(x);
	free(y);
	free(full_x);
	free(matrix);
	return all_pass;
}
int test_dtrsv(unsigned int test_size,unsigned int test_paranoia,double tolerance)
{
	// This function uses dcopy() and dtrmv(). It is
	// assumed that they passed. 
	// the test.
	int all_pass = 1;
	int pass = 1;
	// run tests on all triangular matrix types
	// test unit increment matrix vector product in X
	unsigned int stride = 1;
	double* matrix = (double*)malloc(test_size*test_size*sizeof(double));
	double* x = (double*)malloc(test_size*stride*sizeof(double));
	double* y = (double*)malloc(test_size*stride*sizeof(double));
	double* b = (double*)malloc(test_size*stride*sizeof(double));
	int matrix_types[8] = {1,5,2,6,3,7,4,8};
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		for(unsigned int j = 0 ; j < 8 ; j++)
		{
			if(matrix_types[j] == 1)
			{
				matrix_random_populate(test_size,test_size,matrix,1);
				make_upper_triangular(test_size,matrix);
			}
			else if(matrix_types[j] == 3)
			{
				matrix_random_populate(test_size,test_size,matrix,1);
				make_lower_triangular(test_size,matrix);
			}
			if((matrix_types[j]%2) == 0)			make_unit_diagonal(test_size,matrix);
			vector_random_populate(test_size*stride,x);
			dcopy(test_size,stride,x,stride,b);
			dtrsv(test_size,matrix,stride,x,matrix_types[j]);
			dtrmv(test_size,matrix,stride,x,stride,y,matrix_types[j]);
			pass = compare_vectors(test_size,stride,y,stride,b,1.0,0,0,tolerance);
			if(!pass)			printf("unit stride type %d dtrsv() test failed\n",matrix_types[j]);
			all_pass &= pass;
		}
	}
	free(x);
	free(y);
	free(b);
	// test strided increment matrix vector product in X, generate a random number 
	// for strides but keep it under 10.
	while(stride <= 1)
	{
		stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	x = (double*)malloc(test_size*stride*sizeof(double));
	y = (double*)malloc(test_size*stride*sizeof(double));
	b = (double*)malloc(test_size*stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		for(unsigned int j = 0 ; j < 8 ; j++)
		{
			if(matrix_types[j] == 1)
			{
				matrix_random_populate(test_size,test_size,matrix,1);
				make_upper_triangular(test_size,matrix);
			}
			else if(matrix_types[j] == 3)
			{
				matrix_random_populate(test_size,test_size,matrix,1);
				make_lower_triangular(test_size,matrix);
			}
			if((matrix_types[j]%2) == 0)			make_unit_diagonal(test_size,matrix);
			vector_random_populate(test_size*stride,x);
			dcopy(test_size,stride,x,stride,b);
			dtrsv(test_size,matrix,stride,x,matrix_types[j]);
			dtrmv(test_size,matrix,stride,x,stride,y,matrix_types[j]);
			pass = compare_vectors(test_size,stride,y,stride,b,1.0,0,0,tolerance);
			if(!pass)			printf("strided type %d dtrsv() test failed\n",matrix_types[j]);
			all_pass &= pass;
		}
	}
	free(x);
	free(y);
	free(b);
	free(matrix);
	return all_pass;
}
int test_ddgmv(unsigned int test_size,unsigned int test_paranoia,double tolerance)
{
	// This function uses dgemv(). It is assumed that it passed 
	// the test.
	int all_pass = 1;
	int pass = 1;
	// run tests on all triangular matrix types
	// test unit increment matrix vector product in X
	unsigned int x_stride = 1;
	unsigned int y_stride = 1;
	double* full_A = (double*)malloc(test_size*test_size*sizeof(double));
	double* diagonal_A = (double*)malloc(test_size*sizeof(double));
	double* x = (double*)malloc(test_size*x_stride*sizeof(double));
	double* y = (double*)malloc(test_size*y_stride*sizeof(double));
	double* full_x = (double*)malloc(test_size*x_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		matrix_random_populate(test_size,test_size,full_A,0);
		make_upper_triangular(test_size,full_A);
		make_lower_triangular(test_size,full_A);
		for(unsigned int j = 0 ; j < test_size ; j++)
		{
			diagonal_A[j] = full_A[j*test_size + j];
		}
		vector_random_populate(test_size*x_stride,x);
		vector_random_populate(test_size*y_stride,y);
		set_value(test_size*x_stride,full_x,0.0);
		dgemv(test_size,test_size,1.0,full_A,x_stride,x,1.0,x_stride,full_x,1);
		// test matrix type 1
		ddgmv(test_size,full_A,x_stride,x,y_stride,y,1);
		pass = compare_vectors(test_size,x_stride,full_x,y_stride,y,1.0,0,0,tolerance);
		if(!pass)			printf("unit stride type 1 ddgmv() test failed\n");
		all_pass &= pass;
		// test matrix type 2
		set_value(test_size*y_stride,y,0.0);
		ddgmv(test_size,diagonal_A,x_stride,x,y_stride,y,2);
		pass = compare_vectors(test_size,x_stride,full_x,y_stride,y,1.0,0,0,tolerance);
		if(!pass)			printf("unit stride type 2 ddgmv() test failed\n");
		all_pass &= pass;
	}
	free(x);
	free(y);
	free(full_x);
	// test strided increment matrix vector product in X and Y, generate a random number 
	// for strides but keep it under 10.
	while(x_stride <= 1)
	{
		x_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	while(y_stride <= 1)
	{
		y_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
	}
	x = (double*)malloc(test_size*x_stride*sizeof(double));
	y = (double*)malloc(test_size*y_stride*sizeof(double));
	full_x = (double*)malloc(test_size*x_stride*sizeof(double));
	for(unsigned int i = 0 ; i < test_paranoia ; i++)
	{
		matrix_random_populate(test_size,test_size,full_A,0);
		make_upper_triangular(test_size,full_A);
		make_lower_triangular(test_size,full_A);
		for(unsigned int j = 0 ; j < test_size ; j++)
		{
			diagonal_A[j] = full_A[j*test_size + j];
		}
		vector_random_populate(test_size*x_stride,x);
		vector_random_populate(test_size*y_stride,y);
		set_value(test_size*x_stride,full_x,0.0);
		dgemv(test_size,test_size,1.0,full_A,x_stride,x,1.0,x_stride,full_x,1);
		// test matrix type 1
		ddgmv(test_size,full_A,x_stride,x,y_stride,y,1);
		pass = compare_vectors(test_size,x_stride,full_x,y_stride,y,1.0,0,0,tolerance);
		if(!pass)			printf("strided type 1 ddgmv() test failed\n");
		all_pass &= pass;
		// test matrix type 2
		set_value(test_size*y_stride,y,0.0);
		ddgmv(test_size,diagonal_A,x_stride,x,y_stride,y,2);
		pass = compare_vectors(test_size,x_stride,full_x,y_stride,y,1.0,0,0,tolerance);
		if(!pass)			printf("strided type 2 ddgmv() test failed\n");
		all_pass &= pass;
	}
	free(x);
	free(y);
	free(full_x);
	free(full_A);
	free(diagonal_A);
	return all_pass;
}
void vector_random_populate(unsigned int size,double* x)
{
	double factor = 1.0/(double)RAND_MAX;
	double range = 10.0;
	for(unsigned int i = 0 ; i < size ; i++)
	{
		x[i] = range*(((double)rand())*factor - 0.5);
	}
}
void matrix_random_populate(unsigned int row_count,unsigned int column_count,double* matrix,int non_singular)
{
	double factor = 1.0/(double)RAND_MAX;
	unsigned int size = row_count*column_count;
	double range = 10.0;
	for(unsigned int i = 0 ; i < size ; i++)
	{
		matrix[i] = range*(((double)rand())*factor - 0.5);
	}
	// if the matrix needs to be non-singular, modify its main diagonal elements so 
	// that each of them is larger than the first norm of its row
	if(non_singular)
	{
		double sum = 0.0;
		for(unsigned int i = 0 ; i < row_count ; i++)
		{
			// sum the absolute values of the row entries but skip the main diagonal 
			// element
			sum = dasum(column_count,1,&matrix[i*row_count]);
			sum -= fabs(matrix[i*row_count + i]);
			if(sum > BLAS_tolerance)
			{
				// Make the main diagonal element larger than the first norm of its row.
				// Here, we set the growth factor to a randomly generated real number.
				// The sign of the main diagonal element is maintained.
				double growth_factor = 2.0*(((double)rand())*factor) + 5.0;
				if(matrix[i*row_count + i] < 0.0)			growth_factor = -growth_factor;
				matrix[i*row_count + i] = growth_factor*sum;
			}
			else
			{
				// this is a zero row, make sure that the main diagonal element is not 
				// zero
				while(fabs(matrix[i*row_count + i]) < BLAS_tolerance)
				{
					matrix[i*row_count + i] = range*(((double)rand())*factor - 0.5);
				}
			}
		}
	}
}
void make_upper_triangular(unsigned int size,double* matrix)
{
	for(unsigned int i = 0 ; i < size ; i++)
	{
		for(unsigned int j = 0 ; j < i ; j++)
		{
			matrix[i*size + j] = 0.0;
		}
	}
}
void make_lower_triangular(unsigned int size,double* matrix)
{
	for(unsigned int i = 0 ; i < size ; i++)
	{
		for(unsigned int j = (i + 1) ; j < size ; j++)
		{
			matrix[i*size + j] = 0.0;
		}
	}
}
void make_unit_diagonal(unsigned int size,double* matrix)
{
	double diagonal = 0.0;
	for(unsigned int i = 0 ; i < size ; i++)
	{
		diagonal = matrix[i*size + i];
		for(unsigned int j = 0 ; j < size ; j++)
		{
			matrix[i*size + j] /= diagonal;
		}
	}
}
void set_value(unsigned int size,double* x,double value)
{
	for(unsigned int i = 0 ; i < size ; i++)
	{
		x[i] = value;
	}
}
int compare_vectors(unsigned int size,unsigned int x_inc,const double* x,unsigned int y_inc,const double* y,double y_factor,unsigned int y_shift_inc,const double* y_shift,double tolerance)
{
	int pass = 1;
	if(y_shift != 0)
	{
		for(unsigned int i = 0 ; i < size ; i++)
		{
			if(fabs(x[i*x_inc] - (y[i*y_inc]*y_factor + y_shift[i*y_shift_inc])) > tolerance)
			{
				printf("vector values mismatch : %u : %u : %u : %e : %e : %e\n",size,x_inc,y_inc,x[i*x_inc],y[i*y_inc]*y_factor + y_shift[i*y_shift_inc],x[i*x_inc] - (y[i*y_inc]*y_factor + y_shift[i*y_shift_inc]));
				pass = 0;
			}
		}
	}
	else
	{
		for(unsigned int i = 0 ; i < size ; i++)
		{
			if(fabs(x[i*x_inc] - (y[i*y_inc]*y_factor)) > tolerance)
			{
				printf("vector values mismatch : %u : %u : %u : %e : %e : %e\n",size,x_inc,y_inc,x[i*x_inc],y[i*y_inc]*y_factor,x[i*x_inc]- y[i*y_inc]*y_factor);
				pass = 0;
			}
		}
	}
	return pass;
}


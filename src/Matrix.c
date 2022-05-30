// Matrix.c
// Ahmed M. Hussein (amhussein4@gmail.com)
// 04/27/2022

#include "Matrix.h"
#include "stdlib.h"
#include "string.h"
#include "BLAS.h"
#include "math.h"
#include "Random.h"
#include "stdio.h"

Matrix* create_matrix(unsigned int rows,unsigned int columns)
{
	Matrix* matrix = (Matrix*)malloc(sizeof(Matrix));
	matrix->rows = 0;
	matrix->columns = 0;
	matrix->entries = 0;
	allocate_matrix(matrix,rows,columns);
	return matrix;
}
unsigned int allocate_matrix(Matrix* matrix,unsigned int rows,unsigned int columns)
{
	// reallocate the matrix if it has incorrect dimensions
	unsigned int count = rows*columns;
	if(matrix->rows*matrix->columns != count)
	{
		if(matrix->entries != 0)		free(matrix->entries);
		matrix->entries = (double*)calloc(count,sizeof(double));
	}
	else   					memset(matrix->entries,0,count*sizeof(double));
	matrix->rows = rows;
	matrix->columns = columns;
	return count;
}
void reset_matrix(Matrix* matrix)
{
	if(matrix == 0)					return;
	if(matrix->entries != 0)		free(matrix->entries);
	matrix->entries = 0;
	matrix->rows = 0;
	matrix->columns = 0;
}
void destroy_matrix(Matrix* matrix)
{
	if(matrix == 0)					return;
	reset_matrix(matrix);
	free(matrix);
}
int copy_matrix(const Matrix* source,Matrix* target)
{
	if(source == 0)			return 0;
	if(target == 0)			return 0;
	allocate_matrix(target,source->rows,source->columns);
	dcopy(source->rows*source->columns,1,source->entries,1,target->entries);
	return 1;
}
double mat_get(const Matrix* matrix,unsigned int row,unsigned int column){return matrix->entries[row*matrix->columns + column];}
void mat_set(Matrix* matrix,unsigned int row,unsigned int column,double value){matrix->entries[row*matrix->columns + column] = value;}
void mat_sum(const Matrix* A,const Matrix* B,Matrix* C)
{
	// Compute the sum of matrices A and B and store in C
	// Matrix dimensions must be consistent, no checks are made to ensure that. 
	allocate_matrix(C,A->rows,A->columns);
	unsigned int size = A->rows*A->columns;
	dcopy(size,1,A->entries,1,C->entries);
	daxpy(size,1,B->entries,1.0,1,C->entries);
}
void mat_diff(const Matrix* A,const Matrix* B,Matrix* C)
{
	// Compute the difference of matrices A and B and store in C
	// Matrix dimensions must be consistent, no checks are made to ensure that. 
	allocate_matrix(C,A->rows,A->columns);
	unsigned int size = A->rows*A->columns;
	dcopy(size,1,A->entries,1,C->entries);
	daxpy(size,1,B->entries,-1.0,1,C->entries);
}
void mat_inc(Matrix* A,const Matrix* B)
{
	// Compute the sum of matrices A and B and store in A
	// Matrix dimensions must be consistent, no checks are made to ensure that. 
	daxpy(A->rows*A->columns,1,B->entries,1.0,1,A->entries);
}
void mat_dec(Matrix* A,const Matrix* B)
{
	// Compute the difference of matrices A and B and store in A
	// Matrix dimensions must be consistent, no checks are made to ensure that. 
	daxpy(A->rows*A->columns,1,B->entries,-1.0,1,A->entries);
}
void mat_scale(Matrix* A,double factor)
{
	dscal(A->rows*A->columns,1,A->entries,factor);
}
void mat_mult(const Matrix* A,const Matrix* B,Matrix* C)
{
	// Compute the product of matrices A and B and store the product in C: C = AB
	// Matrix dimensions must be consistent, no checks are made to ensure that. 
	allocate_matrix(C,A->rows,B->columns);
	for(unsigned int i = 0 ; i < B->columns ; i++)
	{
		dgemv(A->rows,A->columns,1.0,A->entries,B->columns,B->entries + i,0.0,B->columns,C->entries + i,1);
	}
}
int mat_solve(const Matrix* A,const Matrix* b,Matrix* x){return mat_solve_jacobi(A,b,x);}
int mat_solve_jacobi(const Matrix* A,const Matrix* b,Matrix* x){return mat_solve_relaxed_jacobi(A,b,x,1.0);}
int mat_solve_relaxed_jacobi(const Matrix* A,const Matrix* b,Matrix* x,double relaxation_factor)
{
	// Solve the system of equations A x = b using the relaxed Jacobi iteration method. The error 
	// is the root mean square error and the relaxation factor is given. The matrix 
	// needs to be diagonally dominant for the iteration to converge but this check 
	// is not made here. 
	unsigned int eq_count = b->rows;
	unsigned int rhs_count = b->columns;
	if(A->rows != A->columns)				return 0;
	if(A->rows != eq_count)					return 0;
	allocate_matrix(x,eq_count,rhs_count);
	Matrix* x_new = create_matrix(eq_count,rhs_count);
	// create matrix y and copy the right hand side to it
	Matrix* y = create_matrix(eq_count,rhs_count);
	// Create the negative inverse diagonal matrix, store the main diagonal only.
	Matrix* inverse_diagonal = create_matrix(eq_count,1);
	for(unsigned int i = 0 ; i < eq_count ; i++)
	{
		inverse_diagonal->entries[i] = -1.0/(A->entries[i*eq_count + i]);
	}
	double tolerance = 1.0e-12;
	double error = 100.0*tolerance;
	int pass = 1;
	unsigned int length = eq_count*rhs_count;
	// The Jacobi iteration: 
	// Decompose matrix A to A = L + D + U where L,D,U are lower traingular, diagonal and 
	// upper triangular matrices respectively
	// x_{i + 1} = D^{-1} (b - (L + U) x_{k}) = D^{-1} (b - (A - D) x_{k})
	// x_{i + 1} = D^{-1} (b - A x_{k}) + D^{-1} D x_{k} = D^{-1} (b - A x_{k}) + x_{k}
	// The relaxed Jacobi iteration:
	// same as Jacobi iteration but x_{i + 1} = x_{k} (1 - f) + (D^{-1} (b - A x_{k}) + x_{k}) f
	// x_{i + 1} = x_{k} + D^{-1} (b - A x_{k})) f
	while(error > tolerance)
	{
		// multiply the current solution x by the matrix A
		mat_mult(A,x,y);
		// subtract the product from the right hand side (or actually subtract the right 
		// hand side from the product since the inverse diagonal matrix is negated)
		daxpy(length,1,b->entries,-1.0,1,y->entries);
		// multiply the negated inverse diagonal matrix by the difference
		if(!dgmat_premult(y,inverse_diagonal,x_new))
		{
			pass = 0;
			break;
		}
		// compute solution increment error
		error = dnrm2(length,1,x_new->entries);
		// add the current solution to solution increment, relax if needed
		if(fabs(relaxation_factor - 1.0) > 1.0e-3)			dscal(length,1,x_new->entries,relaxation_factor);
		daxpy(length,1,x_new->entries,1.0,1,x->entries);
	}
	destroy_matrix(inverse_diagonal);
	destroy_matrix(x_new);
	destroy_matrix(y);
	return pass;
}
int mat_solve_ge(const Matrix* A,const Matrix* b,Matrix* x)
{
	// Solve using Gauss elimination with partial pivoting, very robust but very slow.
	// Matrix dimensions are assumed to be consistent, no checks are made to ensure that. 
	unsigned int eq_count = b->rows;
	unsigned int rhs_count = b->columns;
	if(A->rows != A->columns)				return 0;
	if(A->rows != eq_count)					return 0;
	double* scales = (double*)malloc(eq_count*sizeof(double));
	// get the maximum for each row and place it in the scales vector
	for(unsigned int i = 0 ; i < eq_count ; i++)
	{
		scales[i] = A->entries[i*eq_count + idamax(eq_count,1,A->entries + i*eq_count)];
	}
	// set the working matrix with the row wise scaled version of this matrix
	Matrix* working_matrix = create_matrix(eq_count,eq_count);
	allocate_matrix(x,eq_count,rhs_count);
	for(unsigned int i = 0 ; i < eq_count ; i++)
	{
		daxpy(eq_count,1,A->entries + i*eq_count,1.0/scales[i],1,working_matrix->entries + i*eq_count);
		daxpy(rhs_count,1,b->entries + i*rhs_count,1.0/scales[i],1,x->entries + i*rhs_count);
	}
	free(scales);
	unsigned int swap_index = 0;
	double temp = 0.0;
	for(unsigned int i = 0 ; i < eq_count ; i++)
	{
		swap_index = i + idamax(eq_count - i,eq_count,working_matrix->entries + i*eq_count + i);
		temp = fabs(working_matrix->entries[swap_index*eq_count + i]);
		if(temp < 1.0E-50)
		{
			// matrix is ill-conditioned, cannot solve the system of equations
			destroy_matrix(working_matrix);
			return 0;
		}
		// swap both rows, physically, in memory
		if(swap_index != i)
		{
			dswap(eq_count,1,working_matrix->entries + i*eq_count,1,working_matrix->entries + swap_index*eq_count);
			dswap(rhs_count,1,x->entries + i*rhs_count,1,x->entries + swap_index*rhs_count);
		}
		for(unsigned int j = i + 1 ; j < eq_count; j++)
		{
			temp = -working_matrix->entries[j*eq_count  + i]/working_matrix->entries[i*eq_count + i];
			daxpy(eq_count - i,1,working_matrix->entries + i*eq_count + i,temp,1,working_matrix->entries + j*eq_count + i);
			daxpy(rhs_count,1,x->entries + i*rhs_count,temp,1,x->entries + j*rhs_count);
		}
	}
	// back substitute to get the final solution, do it once per right hand side
	for(unsigned int i = 0 ; i < rhs_count ; i++)
	{	
		dtrsv(eq_count,working_matrix->entries,rhs_count,x->entries + i,1);
	}
	destroy_matrix(working_matrix);
	return 1;
}
int mat_solve_gmres(const Matrix* A,const Matrix* b,Matrix* x)
{
	// Solve the system of equations A x = b using the Generalized Minimum Residual 
	// method by Saad and Schultz. 
	// A is assumed to be a square, non-singular matrix, method works for all types 
	// of non-singular matrices regardless of diagonal dominance, definiteness, sparsity, ...
	// b is a column vector, in case of b having multiple columns, only the first column 
	// will be solved for. To solve for all other columns, separate calls need to be 
	// made to this function
	unsigned int size = A->rows;
	if(size != A->columns)			return 0;
	if(size != b->rows)				return 0;
	// get the norm of the rhs vector
	double b_norm = dnrm(size,b->columns,b->entries);
	// initialize an Arnoldi iteration with the b vector. Arnoldi iterations successively 
	// generate orthonormal column vectors q1, q2, ... and populate an upper Hessenberg 
	// matrix H such that after n iterations, AQn = Qn+1 Hn where Qn is a matrix of the first 
	// n q vectors, Qn+1 is [Qn,qn+1] and Hn is an n+1xn upper Hessenberg matrix
	// Q and H matrices will be created to hold the maximum number of q vectors that will ever 
	// be used in the solution. Solution is guaranteed to converge in less than m iterations 
	// where m is the number of equations and unknowns in the system A x = b
	Matrix* Q = create_matrix(size,size);
	Matrix* H = create_matrix(size,size);	
	Matrix* e = create_matrix(size,1);
	init_arnoldi_iteration(b,Q);
	// the e vector that will be used later is the product of Q^T and b, since the first column 
	// of Q is the normalized b vector and Q is orthonormal, the product Q^T b has only one non-zero 
	// entry in the first component and it is equal to the norm of b. Create and initialize the 
	// product vector e
	mat_set(e,0,0,b_norm);
	// run Arnoldi iteration, every iteration, we get a new q vector and populate a new H column
	// the most we can do is to run m - 1 iterations, there is no need to run the mth iteration
	// Arnoldi iteration order is 0-based with the first call (order = 0) being the initialization 
	// call already made
	double* angle_cosines = (double*)malloc((size - 1)*sizeof(double));
	double* angle_sines = (double*)malloc((size - 1)*sizeof(double));
	double previous_h = 0.0;
	double current_h = 0.0;
	double tolerance = 1.0e-6;
	double error = 0.0;
	// if solution never converges during "size" iterations, then use the maximum convergence_i 
	// which is the maximum number of Givens rotations required to bring H to an upper triangular 
	// form. 
	unsigned int convergence_i = size - 1;
	for(unsigned int i = 0 ; i < size ; i++)
	{
		run_arnoldi_iteration(A,Q,H,i + 1);
		// at this point, Q and H are populated such that AQn = Qn+1 Hn
		// let solution x_n = Qn y where x is mx1 vector and y is nx1 vector
		// the residual vector b - A xn = b - A Qn y = b - Qn+1 Hn y
		// the norm of the residual vector does not change when premultiplied by 
		// the transpose of Qn+1
		// r = ||b - A xn|| = ||b - Qn+1 Hn y|| = ||Qn+1^T b - Hn y||
		// since columns of Q are orthonormal and q1 = b/||b|| (from Arnoldi initialization), 
		// then Qn+1^T b = ||b|| e where e is an n+1 column vector with all zeros except for the 
		// first component
		// r = || ||b|| e - Hn y||, Hn and ||b|| are known at this point
		// Compute Givnes rotation that would eliminate the subdiagonal element of the ith row of the H 
		// matrix and apply it to e only. If applied to both H and e, H will become an upper triangular 
		// matrix and e will be transformed consistently. We do not apply it to H though because 
		// we don't need to do this here as will be shown later. Since Hn is an upper Hessenberg matrix, 
		// we only need to eliminate the zeros along the 1st subdiagonal. 
		// Givens rotation application should not be performed during the last iteration because there 
		// is nothing to eliminate
		if(i < size - 1)
		{
			if(i > 0)
			{
				// for all the columns except for the first one, do the following
				// use the angle cosine and sine from the last Givens rotation (the one for the last column 
				// zero elimination) to transform the H column that has just been generated (in this iteration) 
				// and then use the transformed values (mainly, its main diagonal element) to compute the 
				// new Givens rotation angle cosine and sine.
				// Why ? Because we don't explicitly transform H as we go. In addition, H is generated 
				// incrementally. This means that in column k reduction, we transform rows k and k + 1, which 
				// in principle should change elements (k,k), (k,k + 1), (k,k + 2), ... (k,m) and (k + 1,k), 
				// (k + 1,k + 1), (k + 1,k + 2), ... (k + 1,m) where m is the matrix size. All elements in rows 
				// k and k + 1 and columns j < k are zeros and will not be affected by this rotation. 
				// If we rotate the current H as is, all the columns k + 1 and above do not exist yet and 
				// they will be zeros. The only elements which will be properly rotated are (k,k) and (k+1,k). 
				// While we don't care, for now, about all the elements in rows k and k + 1 and columns k + 2 
				// and above, we particularly care about the element (k + 1,k + 1) because this is the element 
				// which will be used for the rotation angle calculation when reducing the next column. Having 
				// an incorrect column value at this point results in an incorrect rotation for the next column 
				// and the error propagates forward. We don't need to worry about the element (k + 2,k + 1) because 
				// in the rotation of k and k + 1 rows, k + 2 is untouched. 
				// So now we are in column k + 1 rotation, we need to transform element (k + 1,k + 1) based on the 
				// previous rotation angle before we apply this column's (k + 1) rotation. Do this now. If this is 
				// the first iteration, we just go ahead and rotate with what we have in H. No risk here since 
				// it is the first rotation. 
				// The rotation application here is taken from the BLAS drot function
				// The main diagonal element we need to transform is at i,i (which is k + 1,k + 1 in the above comment) 
				// and the other row element (k,k + 1 in the above comment) is at i - 1,i. Notice that this element 
				// at i - 1,i has just been generated in this iteration
				// To add to the complexity, we need to apply all the past rotations to elements of the i column. 
				// Why ? Because the first rotation changes elements 1,j for j >= 1, which includes this column k, 
				// the second rotation changes elements 2,j for j >= 2, which, in general, includes this column, 
				// the third rotation changes elements 3,j for j >= 3, ... etc. This means that in order to get 
				// the correct element k,k, we need to apply the last rotation between k,k and k - 1,k, but in 
				// order to get correct element k - 1,k, we need to apply the second to last rotation between 
				// k - 1,k and k - 2,k, and so on. The angle cosines and sines are stored, use them to do this here
				// previous_h will always be the kth column element in the top row we are transforming
				previous_h = H->entries[i];
				for(unsigned int j = 0 ; j < i ; j++)
				{
					current_h = H->entries[(j + 1)*size + i];
					previous_h = angle_cosines[j]*current_h - angle_sines[j]*previous_h;
				}
				// after all rotations, do not forget to compute Givens rotation based on the transformed values. In 
				// all cases, do not actually apply the rotation on H, do it on the fly and apply the rotations 
				// outside of this loop after convergence. 
				drotg(previous_h,H->entries[(i + 1)*size + i],angle_cosines + i,angle_sines + i);
			}
			else 			drotg(H->entries[i*size + i],H->entries[(i + 1)*size + i],angle_cosines,angle_sines);
			drot(1,1,e->entries + i,1,e->entries + i + 1,angle_cosines[i],angle_sines[i]);
			// now we have, in theory, an upper triangular matrix above a zero row multiplied by y. 
			// The zero row, when multiplied by y vector, gives 0, which, when subtracted from the
			// transformed e column vector in the residual norm to minimize, gives a number independent of y. 
			// Varying y changes the residual norm but this last component acts like a fixed norm component 
			// that cannot be minimized. We need to find y that minimizes the rest of the system. 
			// This is straightforward now since the upper nxn part of the transformed Hn is an 
			// upper triangular matrix. 
			// we will not do this here however, we will only compute the residual norm and if it is below 
			// the acceptable tolerance, we will apply all the Givens rotation on H and e in a separate loop. 
			// This is why we didn't apply Givens rotations on H earlier. 
			// We applied Givens rotations on e however because the last component of e is the residual norm 
			// since after we solve for the upper nxn triangular matrix, all components vanish and the error 
			// becomes the magnitude of the last component
			// Again, this cannot be applied during the last iteration. In general, if we reach the 
			// last iteration, there is no point in computing the error since we will have to apply 
			// Givens rotations all the way anyway. The purpose of error calculations is to break early 
			// if possible
			error = fabs(mat_get(e,i + 1,0));
			// if the error is less than the error tolerance, break out of the loop
			if(error < tolerance)
			{
				convergence_i = i;
				break;
			}
		}
	}
	// Now that we know the iteration at which convergence takes place, apply all 
	// givens rotations up to this iteration then solve the system of equations. Givens 
	// rotations need to be applied on both the H matrix and the e vector. Start with the 
	// original, unrotated, e vector. This transforms H to an upper triangular matrix. We 
	// only need to apply 1 rotation per column since H is an upper Hessenberg matrix. 
	//convergence_i = 2;
	// Solution size is the number of rows/column of the upper triangular matrix resulting from 
	// H matrix rotation. It is 1 more than the convergence_i because this is where we stopped 
	// rotating H (in principle)
	unsigned int solution_size = convergence_i + 1;
	mat_zero(e);
	mat_set(e,0,0,b_norm);
	// make sure we don't do a rotation for the last column in all cases
	unsigned int rotation_count = ((solution_size < size) ? solution_size : size - 1);
	for(unsigned int i = 0 ; i < rotation_count ; i++)
	{
		// there is no need to recompute the angle sines and cosines, the correct ones were calculated 
		// during the iteration above
		drot(size,1,H->entries + i*size,1,H->entries + (i + 1)*size,angle_cosines[i],angle_sines[i]);
		drot(1,1,e->entries + i,1,e->entries + i + 1,angle_cosines[i],angle_sines[i]);
	}
	// solve the kxk upper triangular matrix H where k is convergence_i
	// before we do that, we need to reshape H to fit kxk rather than size x size
	Matrix* T = create_matrix(solution_size,solution_size);
	mat_reshape(H,T);
	dtrsv(solution_size,T->entries,1,e->entries,1);
	// form x = Q y where y is stored in e. Only use enough columns of Q in the product and no more. Q is 
	// always 1 column more than needed. Regular matrix multiplication doesn't work here and is generally 
	// inefficient
	allocate_matrix(x,size,1);
	for(unsigned int i = 0 ; i < solution_size ; i++)
	{
		daxpy(size,size,Q->entries + i,e->entries[i],1,x->entries);
	}
	destroy_matrix(Q);
	destroy_matrix(H);
	destroy_matrix(T);
	destroy_matrix(e);
	free(angle_cosines);
	free(angle_sines);
	return 1;
}
int dgmat_premult(const Matrix* matrix,const Matrix* diagonal,Matrix* product)
{
	// the passed diagonal matrix is assumed to be given as an nx1 vector not an nxn diagonal matrix
	// check that the matrix is a square matrix and that the dimensions are consistent
	if(diagonal->rows != matrix->rows)		return 0;
	allocate_matrix(product,matrix->rows,matrix->columns);
	for(unsigned int i = 0 ; i < matrix->columns ; i++)
	{
		ddgmv(diagonal->rows,diagonal->entries,matrix->columns,matrix->entries + i,matrix->columns,product->entries + i,2);
	}
	return 1;
}
int dgmat_postmult(const Matrix* matrix,const Matrix* diagonal,Matrix* product)
{
	// the passed diagonal matrix is assumed to be given as an nx1 vector not an nxn diagonal matrix
	// check that the matrix is a square matrix and that the dimensions are consistent
	if(matrix->columns != diagonal->rows)		return 0;
	allocate_matrix(product,matrix->rows,matrix->columns);
	for(unsigned int i = 0 ; i < matrix->columns ; i++)
	{
		daxpy(matrix->rows,matrix->columns,matrix->entries + i,diagonal->entries[i],matrix->columns,product->entries + i);
	}
	return 1;
}
int mat_inv(const Matrix* matrix,Matrix* inverse)
{
	// check that the matrix is a square matrix
	if(matrix->rows != matrix->columns)			return 0;
	// create RHS matrix
	Matrix* rhs = create_matrix(matrix->rows,matrix->columns);
	for(unsigned int i = 0 ; i < matrix->rows ; i++)
	{
		rhs->entries[i*matrix->columns + i] = 1.0;
	}
	int pass = mat_solve_ge(matrix,rhs,inverse);
	destroy_matrix(rhs);
	if(!pass)				destroy_matrix(inverse);
	return pass;
}
void mat_rand(Matrix* matrix,double start,double end,int diagonally_dominant)
{
	// populates the passed matrix with random values given its dimensions
	unsigned int size = matrix->rows*matrix->columns;
	for(unsigned int i = 0 ; i < size ; i++)
	{
		matrix->entries[i] = rand_uniform_interval(start,end);
	}
	// guarantee diagonal dominance if needed
	if(diagonally_dominant)
	{
		if(matrix->rows != matrix->columns)					return;
		if(matrix->rows == 1)								return;
		double sum = 0.0;
		for(unsigned int i = 0 ; i < matrix->rows ; i++)
		{
			// sum all absolute values of row entries
			sum = dasum(matrix->columns,1,matrix->entries + i*matrix->columns);
			// subtract the main diagonal entry
			sum -= matrix->entries[i*matrix->columns + i];
			// update the main diagonal term by scaling the sum with a factor that 
			// is greater than one
			matrix->entries[i*matrix->columns + i] = rand_uniform_interval(1.0,10.0)*sum;
		}
	}
}
void mat_nrand(Matrix* matrix,double mean,double standard_deviation)
{
	// populates the passed matrix with normally distributed random values given its dimensions
	unsigned int size = matrix->rows*matrix->columns;
	for(unsigned int i = 0 ; i < size ; i++)
	{
		matrix->entries[i] = rand_normal(mean,standard_deviation);
	}
}
double mat_norm2(const Matrix* matrix)
{
	return dnrm2(matrix->rows*matrix->columns,1,matrix->entries);
}
void mat_print(const Matrix* matrix)
{
	for(unsigned int i = 0 ; i < matrix->rows ; i++)
	{
		for(unsigned int j = 0 ; j < matrix->columns ; j++)
		{
			printf("%e\t",matrix->entries[i*matrix->columns + j]);
		}
		printf("\n");
	}
}
void mat_zero(Matrix* matrix)
{
	memset(matrix->entries,0,matrix->rows*matrix->columns*sizeof(double));
}
void mat_identity(Matrix* matrix)
{
	if(matrix->rows != matrix->columns)		return;
	mat_zero(matrix);
	for(unsigned int i = 0 ; i < matrix->rows ; i++)
	{
		matrix->entries[i*matrix->columns + i] = 1.0;
	}
}
void mat_transpose(const Matrix* matrix,Matrix* transpose)
{
	allocate_matrix(transpose,matrix->columns,matrix->rows);
	for(unsigned int i = 0 ; i < matrix->rows ; i++)
	{
		dcopy(matrix->columns,1,matrix->entries + i*matrix->columns,matrix->rows,transpose->entries + i);
	}
}
int init_arnoldi_iteration(const Matrix* b,Matrix* Q)
{
	// This function initializes Arnoldi iteration using a column vector b
	// Matrix Q is the matrix of orthonormal vectors, it is assumed that this matrix 
	// has been properly allocated before calling this function
	// b can have any number of columns, only the first column will be used for 
	// initialization
	unsigned int rows = b->rows;
	if(Q->rows != rows)			return 0;
	double norm = dnrm(rows,b->columns,b->entries);
	dcopy(rows,b->columns,b->entries,Q->columns,Q->entries);
	dscal(rows,Q->columns,Q->entries,1.0/norm);
	return 1;
}
int run_arnoldi_iteration(const Matrix* A,Matrix* Q,Matrix* H,unsigned int order)
{
	// This function performs Arnoldi iteration of the given order based on matrix A
	// The matrices Q and H are the orthogonormal vectors matrix and the Hessenberg matrix 
	// respectively. They are assumed to have been properly allocated before calling this function
	// The order is the order of the orthonormal vector that will be generated at the 
	// end of the call to this function
	// The dimension of the orthonormal vectors will be taken from the A matrix and 
	// Q matrix needs to have enough space to store the new vector generated in 
	// this iteration
	// Order is 0-based, an order 0 call is an initialization call which should have 
	// been made before calling this function by calling init_arnoldi_iteration()
	if(order == 0)				return 0;
	unsigned int rows = A->rows;
	unsigned int columns = A->columns;
	unsigned int stride = Q->columns;
	// A needs to be a square matrix
	if(rows != columns)			return 0;
	if(Q->rows != columns)		return 0;
	if(stride < order)			return 0;
	if(order == rows)
	{
		// if this is the last iteration, then the Q matrix is already full, we only need to 
		// populate the last column of the H vector
		// Since A = Q H Qt --> H = Qt A Q, Q is known and all columns of H, except for the last 
		// one, are known
		// last H column = Qt A qm where qm is the last column of Q
		// <ultiply A by qm and store in a temporary column vector. We cannot use the last 
		// column of H because it will be overwritten in the next multiplication step
		Matrix* v = create_matrix(rows,1);
		dgemv(rows,columns,1.0,A->entries,stride,Q->entries + order - 1,0.0,1,v->entries,1);
		// premultiply the result by Qt 
		dgemv(rows,columns,1.0,Q->entries,1,v->entries,0.0,stride,H->entries + order - 1,2);
		destroy_matrix(v);
		return 1;
	}
	double factor = 0.0;
	// (order - 1) and (order) in the call below because column indices are zero-based and order is zero-based
	dgemv(rows,columns,1.0,A->entries,stride,Q->entries + order - 1,0.0,stride,Q->entries + order,1);
	for(unsigned int i = 0 ; i < order ; i++)
	{
		factor = ddot(rows,stride,Q->entries + i,stride,Q->entries + order);
		daxpy(rows,stride,Q->entries + i,-factor,stride,Q->entries + order);
		mat_set(H,i,order - 1,factor);
	}
	factor = dnrm(rows,stride,Q->entries + order);
	dscal(rows,stride,Q->entries + order,1.0/factor);
	mat_set(H,order,order - 1,factor);
	return 1;
}
void mat_reshape(const Matrix* source,Matrix* target)
{
	// This function reshapes the entries of source matrix to fit in 
	// target matrix. The sizes of the source and target matrices indicate 
	// the reshaping procedure so there is no need to pass matrix sizes 
	// independently. All entries in the target matrix are overwritten and 
	// all extra entries in the source matrix are discarded. 
	// Rows of source matrix are copied to rows of target matrix. If the 
	// width of source is more than width of target, the remaining part of 
	// each source row is discarded, otherwise, the remaining part of each 
	// target row is padded with zeros. This happens to each source row, if 
	// the number of source rows is more than the number of target rows, the 
	// ramining rows are discarded, otherwise, the remaining target rows are 
	// padded with zeros
	mat_zero(target);
	unsigned int height = ((target->rows < source->rows) ? target->rows : source->rows);
	unsigned int width = ((target->columns < source->columns) ? target->columns : source->columns);
	for(unsigned int i = 0 ; i < height ; i++)
	{
		dcopy(width,1,source->entries + i*source->columns,1,target->entries + i*target->columns);
	}
}
int test_mat(unsigned int size,double tolerance,unsigned int paranoia)
{
	// use a random number of right hand sides
	unsigned int rhs_count = rand_int(1,10);
	// create all required matrices
	Matrix* A = create_matrix(size,size);
	Matrix* b = create_matrix(size,rhs_count);
	Matrix* x = create_matrix(size,rhs_count);
	Matrix* C = create_matrix(size,rhs_count);
	Matrix* Ainv = create_matrix(size,size);
	Matrix* D = create_matrix(size,size);
	Matrix* I = create_matrix(size,size);
	mat_identity(I);
	int pass = 0;
	for(unsigned int i = 0 ; i < paranoia ; i++)
	{
		// generate system of equations
		mat_rand(A,-10.0,10.0,1);
		mat_rand(b,-5.0,5.0,0);
		// test general solution
		pass = mat_solve(A,b,x);
		if(pass)
		{
			mat_mult(A,x,C);
			mat_dec(C,b);
			printf("solution error = %e\n",mat_norm2(C));
		}
		else		printf("solution of linear system of equations of size %u failed\n",size);
		// test GE solution
		pass = mat_solve_ge(A,b,x);
		if(pass)
		{
			
			mat_mult(A,x,C);
			mat_dec(C,b);
			printf("GE solution error = %e\n",mat_norm2(C));
		}
		else		printf("GE solution of linear system of equations of size %u failed\n",size);
		// test matrix inversion
		pass = mat_inv(A,Ainv);
		if(pass)
		{
			mat_mult(A,Ainv,D);
			mat_dec(D,I);
			printf("inversion error = %e\n",mat_norm2(D));
		}
		else		printf("matrix inversion failed\n");
	}
	destroy_matrix(A);
	destroy_matrix(b);
	destroy_matrix(x);
	destroy_matrix(C);
	destroy_matrix(Ainv);
	destroy_matrix(D);
	destroy_matrix(I);
	return pass;
}

/*namespace EZ
{
	namespace Math
	{
		
		void Matrix::SetRow(const unsigned int& row_index,double* row_entries)
		{
			// Set the row of index row_index with the entries in row_entries
			// The array row_entries is assumed to have been properly allocated 
			// and populated. No further checks are made here. 
			BLAS::Copy(column_count,1,row_entries,1,&entries[row_index*column_count]);
		}
		double Matrix::operator()(const unsigned int& row_index) const{return entries[row_index*column_count];}
		double Matrix::operator()(const unsigned int& row_index,const unsigned int& column_index) const{return entries[row_index*column_count + column_index];}
		void Matrix::operator()(const unsigned int& row_index,const unsigned int& column_index,const double& value){entries[row_index*column_count + column_index] = value;}
		void Matrix::Increment(const unsigned int& row_index,const unsigned int& column_index,const double& value){entries[row_index*column_count + column_index] += value;}
		void Matrix::Decrement(const unsigned int& row_index,const unsigned int& column_index,const double& value){entries[row_index*column_count + column_index] -= value;}


		Matrix Matrix::operator*(const double& factor) const
		{
			Matrix product(*this);
			BLAS::Scale(row_count*column_count,1,product.entries,factor);
			return product;
		}
		double Matrix::operator^(const Matrix& matrix) const
		{
			// this is a contraction operator, both matrices should have 
			// the same number of rows and columns
			if(matrix.row_count != row_count)			return 0.0;
			if(matrix.column_count != column_count)		return 0.0;
			double contraction = 0.0;
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				for(unsigned int j = 0 ; j < column_count ; j++)
				{
					contraction += (entries[i*column_count + j]*matrix.entries[i*column_count + j]);
				}
			}
			return contraction;
		}
		void Matrix::SumRows(double* sums) const
		{
			double sum = 0.0;
			unsigned int index = 0;
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				sum = 0.0;
				for(unsigned int j = 0 ; j < column_count ; j++)
				{
					sum += entries[index++];
				}
				sums[i] = sum;
			}
		}
		bool Matrix::WriteRow(const unsigned int& row,FILE* file) const
		{
			return (column_count == fwrite(&entries[row*column_count],sizeof(double),column_count,file));
		}
		bool Matrix::ReadRow(const unsigned int& row,FILE* file) const
		{
			return (column_count == fread(&entries[row*column_count],sizeof(double),column_count,file));
		}
		const double* Matrix::Row(const unsigned int& row_index) const{return &entries[row_index*column_count];}
		const double* Matrix::Entries() const{return entries;}
		void Matrix::Initialize(){entries = 0;}

		SparseMatrix::SparseMatrix(){Initialize();}
		SparseMatrix::SparseMatrix(const SparseMatrix& matrix) : BaseMatrix(matrix)
		{
			*this = matrix;
		}
		SparseMatrix::SparseMatrix(const unsigned int& size)
		{
			Initialize();
			Allocate(size,size);
		}
		SparseMatrix::SparseMatrix(const unsigned int& target_row_count,const unsigned int& target_column_count)
		{
			Initialize();
			Allocate(target_row_count,target_column_count);
		}
		SparseMatrix::~SparseMatrix(){Reset();}
		SparseMatrix& SparseMatrix::operator=(const SparseMatrix& matrix)
		{
			Allocate(matrix.row_count,matrix.column_count);
			for(std::map<unsigned int,double>::const_iterator entry = matrix.entries.begin() ; entry != matrix.entries.end() ; entry++)
			{
				entries[entry->first] = entry->second;
			}
			return *this;
		}
		void SparseMatrix::Reset()
		{
			entries.clear();
			BaseMatrix::Reset();
			Initialize();
		}
		void SparseMatrix::Allocate(const unsigned int& target_row_count,const unsigned int& target_column_count)
		{
			Reset();
			row_count = target_row_count;
			column_count = target_column_count;
			entries.clear();
		}
		void SparseMatrix::Randomize()
		{
			entries.clear();
			double density = 0.01;
			unsigned int size = row_count*column_count;
			unsigned int entry_count = (unsigned int)floor(density*size + 0.5);
			for(unsigned int i = 0 ; i < entry_count ; i++)
			{
				entries[Random::UniformInteger(0,size)] = Random::Uniform(-10.0,10.0);
			}
		}
		void SparseMatrix::RandomizeDiagonallyDominant()
		{
			Randomize();
			double* row_sums = new double[row_count];
			double* diagonals = new double[row_count];
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				row_sums[i] = 0.0;
				diagonals[i] = 0.0;
			}
			// sum the entries per row
			unsigned int row = 0;
			for(std::map<unsigned int,double>::const_iterator entry = entries.begin() ; entry != entries.end() ; entry++)
			{
				row = entry->first/column_count;
				row_sums[row] += fabs(entry->second);
				if(row == (entry->first%column_count))		diagonals[row] = entry->second;
			}
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				// subtract the main diagonal entry from the row sum
				row_sums[i] -= diagonals[i];
				// update the main diagonal term by scaling the sum with a factor that 
				// is greater than one
				if(row_sums[i] < 1.0e-3)			row_sums[i] = 1.0;
				entries[i*column_count + i] = Random::Uniform(1.0,10.0)*row_sums[i];
			}
		}
		double SparseMatrix::operator()(const unsigned int& row_index,const unsigned int& column_index) const
		{
			unsigned int index = row_index*column_count + column_index;
			std::map<unsigned int,double>::const_iterator entry = entries.find(index);
			if(entry == entries.end())		return 0.0;
			return entry->second;
		}
		void SparseMatrix::operator()(const unsigned int& row_index,const unsigned int& column_index,const double& value)
		{
			unsigned int index = row_index*column_count + column_index;
			entries[index] = value;
		}
		Matrix SparseMatrix::operator*(const Matrix& matrix) const
		{
			unsigned int row = 0;
			unsigned int column = 0;
			unsigned int matrix_columns = matrix.ColumnCount();
			Matrix product(row_count,matrix_columns);
			for(std::map<unsigned int,double>::const_iterator entry = entries.begin() ; entry != entries.end() ; entry++)
			{
				row = entry->first/column_count;
				column = entry->first%column_count;
				for(unsigned int j = 0 ; j < matrix_columns ; j++)
				{
					product.Increment(row,j,entry->second*matrix(column,j));
				}
			}
			return product;
		}
		void SparseMatrix::SumRows(double* sums) const
		{
			for(std::map<unsigned int,double>::const_iterator entry = entries.begin() ; entry != entries.end() ; entry++)
			{
				sums[entry->first/column_count] += entry->second;
			}
		}
		void SparseMatrix::Initialize(){entries.clear();}
	}
}*/


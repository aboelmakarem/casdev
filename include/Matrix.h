// Matrix.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 04/27/2022

#ifndef MATRIX_H_
#define MATRIX_H_

typedef struct Matrix
{
	unsigned int rows;
	unsigned int columns;
	double* entries;
} Matrix;

Matrix* create_matrix(unsigned int rows,unsigned int columns);
unsigned int allocate_matrix(Matrix* matrix,unsigned int rows,unsigned int columns);
void reset_matrix(Matrix* matrix);
void destroy_matrix(Matrix* matrix);
int copy_matrix(const Matrix* source,Matrix* target);
double mat_get(const Matrix* matrix,unsigned int row,unsigned int column);
void mat_set(Matrix* matrix,unsigned int row,unsigned int column,double value);
void mat_sum(const Matrix* A,const Matrix* B,Matrix* C);
void mat_diff(const Matrix* A,const Matrix* B,Matrix* C);
void mat_inc(Matrix* A,const Matrix* B);
void mat_dec(Matrix* A,const Matrix* B);
void mat_scale(Matrix* A,double factor);
void mat_mult(const Matrix* A,const Matrix* B,Matrix* C);
int mat_solve(const Matrix* A,const Matrix* b,Matrix* x);
int mat_solve_jacobi(const Matrix* A,const Matrix* b,Matrix* x);
int mat_solve_relaxed_jacobi(const Matrix* A,const Matrix* b,Matrix* x,double relaxation_factor);
int mat_solve_ge(const Matrix* A,const Matrix* b,Matrix* x);
int mat_solve_gmres(const Matrix* A,const Matrix* b,Matrix* x);
int dgmat_premult(const Matrix* matrix,const Matrix* diagonal,Matrix* product);
int dgmat_postmult(const Matrix* matrix,const Matrix* diagonal,Matrix* product);
int mat_inv(const Matrix* matrix,Matrix* inverse);
void mat_rand(Matrix* matrix,double start,double end,int diagonally_dominant);
void mat_nrand(Matrix* matrix,double mean,double standard_deviation);
double mat_norm2(const Matrix* matrix);
void mat_print(const Matrix* matrix);
void mat_zero(Matrix* matrix);
void mat_identity(Matrix* matrix);
void mat_transpose(const Matrix* matrix,Matrix* transpose);
int init_arnoldi_iteration(const Matrix* b,Matrix* Q);
int run_arnoldi_iteration(const Matrix* A,Matrix* Q,Matrix* H,unsigned int order);
void mat_reshape(const Matrix* source,Matrix* target);
int test_mat(unsigned int size,double tolerance,unsigned int paranoia);

/*namespace EZ
{
	namespace Math
	{
		class Matrix : public BaseMatrix
		{
		public:
			void SetRow(const unsigned int& row_index,double* row_entries);
			double operator()(const unsigned int& row_index) const;
			double operator()(const unsigned int& row_index,const unsigned int& column_index) const;
			void operator()(const unsigned int& row_index,const unsigned int& column_index,const double& value);
			void Increment(const unsigned int& row_index,const unsigned int& column_index,const double& value);
			void Decrement(const unsigned int& row_index,const unsigned int& column_index,const double& value);
			Matrix operator*(const double& factor) const;
			double operator^(const Matrix& matrix) const;
			void SumRows(double* sums) const;
			const double* Entries() const;
		};

		class SparseMatrix : public BaseMatrix
		{
		public:
			SparseMatrix();
			SparseMatrix(const SparseMatrix& matrix);
			SparseMatrix(const unsigned int& size);
			SparseMatrix(const unsigned int& target_row_count,const unsigned int& target_column_count);
			~SparseMatrix();
			SparseMatrix& operator=(const SparseMatrix& matrix);
			void Reset();
			void Allocate(const unsigned int& target_row_count,const unsigned int& target_column_count);
			void Randomize();
			void RandomizeDiagonallyDominant();
			double operator()(const unsigned int& row_index,const unsigned int& column_index) const;
			void operator()(const unsigned int& row_index,const unsigned int& column_index,const double& value);
			Matrix operator*(const Matrix& matrix) const;
			void SumRows(double* sums) const;
			
		private:
			void Initialize();
			std::map<unsigned int,double> entries;
		};
	}
}*/

#endif


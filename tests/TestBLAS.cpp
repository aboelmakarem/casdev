/*	TestBLAS.cpp
	Ahmed M. Hussein (amhussein4@gmail.com)
	05/01/2022

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
}

int main(int argc,char** argv)
{
	if(argc < 3)
	{
		printf("error: missing run arguments\n");
		printf("usage:testblas test_size test_count\n");
		return 1;
	}
	unsigned int test_size = atoi(argv[1]);
	unsigned int test_count = atoi(argv[2]);
	printf("running %u BLAS tests of size %u\n",test_count,test_size);
	if(test_BLAS(test_size,test_count,1.0e-8,1))		printf("All BLAS tests passed\n");
	else 												printf("Some or all BLAS tests failed\n");
	return 0;
}


/*	Complex.c

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

#include "Complex.h"
#include "stdlib.h"
#include "math.h"

Complex* create_complex()
{
	Complex* complex = (Complex*)malloc(sizeof(Complex));
	complex->real = 0.0;
	complex->imaginary = 0.0;
	return complex;
}
void reset_complex(Complex* complex)
{
	complex->real = 0.0;
	complex->imaginary = 0.0;
}
void destroy_complex(Complex* complex)
{
	if(complex == 0)		return;
	reset_complex(complex);
	free(complex);
}
void copy_complex(const Complex* source,Complex* target)
{
	target->real = source->real;
	target->imaginary = source->imaginary;
}
void complex_amp_phase(Complex* complex,double amplitude,double phase)
{
	complex->real = amplitude*cos(phase);
	complex->imaginary = amplitude*sin(phase);
}
double complex_amp(const Complex* complex){return sqrt(complex_amp2(complex));}
double complex_phase(const Complex* complex){return atan2(complex->imaginary,complex->real);}
double complex_amp2(const Complex* complex){return (complex->real*complex->real + complex->imaginary*complex->imaginary);}
int complex_is_real(const Complex* complex)
{
	if(fabs(complex->imaginary/complex->real) < 1.0e-12)		return 1;
	return 0;
}
int complex_is_imaginary(const Complex* complex)
{
	if(fabs(complex->real/complex->imaginary) < 1.0e-12)		return 1;
	return 0;
}
void complex_conjugate(const Complex* complex,Complex* conjugate)
{
	conjugate->real = complex->real;
	conjugate->imaginary = -complex->imaginary;
}
void complex_add(const Complex* a,const Complex* b,Complex* c)
{
	c->real = a->real + b->real;
	c->imaginary = a->imaginary + b->imaginary;
}
void complex_inc(Complex* a,const Complex* b)
{
	a->real += b->real;
	a->imaginary += b->imaginary;
}
void complex_diff(const Complex* a,const Complex* b,Complex* c)
{
	c->real = a->real - b->real;
	c->imaginary = a->imaginary - b->imaginary;
}
void complex_dec(Complex* a,const Complex* b)
{
	a->real -= b->real;
	a->imaginary -= b->imaginary;
}
void complex_mult(const Complex* a,const Complex* b,Complex* c)
{
	c->real = a->real*b->real - a->imaginary*b->imaginary;
	c->imaginary = a->real*b->imaginary + a->imaginary*b->real;
}
void complex_divide(const Complex* a,const Complex* b,Complex* c)
{
	Complex* bconj = create_complex();
	Complex* prod = create_complex();
	complex_conjugate(b,bconj);
	complex_mult(a,bconj,prod);
	complex_scale(prod,1.0/complex_amp2(b),c);
	destroy_complex(bconj);
	destroy_complex(prod);
}
void complex_scale(const Complex* a,double factor,Complex* c)
{
	c->real = a->real*factor;
	c->imaginary = a->imaginary*factor;
}


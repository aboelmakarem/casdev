/*	Complex.h

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

#ifndef COMPLEX_H_
#define COMPLEX_H_

typedef struct Complex
{
	double real;
	double imaginary;
} Complex;

Complex* create_complex();
void reset_complex(Complex* complex);
void destroy_complex(Complex* complex);
void copy_complex();
void complex_amp_phase(Complex* complex,double amplitude,double phase);
double complex_amp(const Complex* complex);
double complex_phase(const Complex* complex);
double complex_amp2(const Complex* complex);
int complex_is_real(const Complex* complex);
int complex_is_imaginary(const Complex* complex);
void complex_conjugate(const Complex* complex,Complex* conjugate);
void complex_add(const Complex* a,const Complex* b,Complex* c);
void complex_inc(Complex* a,const Complex* b);
void complex_diff(const Complex* a,const Complex* b,Complex* c);
void complex_dec(Complex* a,const Complex* b);
void complex_mult(const Complex* a,const Complex* b,Complex* c);
void complex_divide(const Complex* a,const Complex* b,Complex* c);
void complex_scale(const Complex* a,double factor,Complex* c);

#endif


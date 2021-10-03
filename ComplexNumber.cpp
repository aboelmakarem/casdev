// ComplexNumber.cpp
// Copyright (c) 2010 Ahmed M. Hussein (amhussein4@gmail.com)
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "ComplexNumber.h"
#include "math.h"

namespace CASDev
{
	ComplexNumber::ComplexNumber(){Initialize();}
	ComplexNumber::ComplexNumber(const ComplexNumber& oNumber){*this = oNumber;}
	ComplexNumber::ComplexNumber(const double& real_part,const double& imaginary_part){Set(real_part,imaginary_part);}
	ComplexNumber::~ComplexNumber(){Reset();}
	ComplexNumber& ComplexNumber::operator=(const ComplexNumber& number)
	{
		real = number.real;
		imaginary = number.imaginary;
		return *this;
	}
	void ComplexNumber::Reset(){Initialize();}
	void ComplexNumber::Set(const double& real_part,const double& imaginary_part)
	{
		real = real_part;
		imaginary = imaginary_part;
	}
	void ComplexNumber::AmplitudePhase(const double& amplitude,const double& phase)
	{
		real = amplitude*cos(phase);
		imaginary = amplitude*sin(phase);
	}
	void ComplexNumber::Real(const double& value){real = value;}
	void ComplexNumber::Imaginary(const double& value){imaginary = value;}
	double ComplexNumber::Real() const{return real;}
	double ComplexNumber::Imaginary() const{return imaginary;}
	double ComplexNumber::Amplitude() const{return sqrt(SquaredAmplitude());}
	double ComplexNumber::Phase() const{return atan2(imaginary,real);}
	double ComplexNumber::SquaredAmplitude() const{return (real*real + imaginary*imaginary);}
	bool ComplexNumber::IsReal(const double& tolerance) const{return (fabs(imaginary) < tolerance);}
	bool ComplexNumber::IsImaginary(const double& tolerance) const{return (fabs(real) < tolerance);}
	ComplexNumber ComplexNumber::operator+(const ComplexNumber& number) const
	{
		return ComplexNumber(real + number.real,imaginary + number.imaginary);
	}
	ComplexNumber ComplexNumber::operator+=(const ComplexNumber& number)
	{
		real += number.real;
		imaginary += number.imaginary;
		return *this;
	}
	ComplexNumber ComplexNumber::operator-(const ComplexNumber& number) const
	{
		return ComplexNumber(real - number.real,imaginary - number.imaginary);
	}
	ComplexNumber ComplexNumber::operator-=(const ComplexNumber& number)
	{
		real -= number.real;
		imaginary -= number.imaginary;
		return *this;
	}
	ComplexNumber ComplexNumber::operator+(const double& number) const{return ComplexNumber(real + number,imaginary);}
	ComplexNumber ComplexNumber::operator+=(const double& number)
	{
		real += number;
		return *this;
	}
	ComplexNumber ComplexNumber::operator-(const double& number) const{return ComplexNumber(real - number,imaginary);}
	ComplexNumber ComplexNumber::operator-=(const double& number)
	{
		real -= number;
		return *this;
	}
	ComplexNumber ComplexNumber::operator*(const ComplexNumber& number) const
	{
		ComplexNumber result;
		result.Real(real*number.real - imaginary*number.imaginary);
		result.Imaginary(real*number.imaginary + imaginary*number.real);
		return result;
	}
	ComplexNumber ComplexNumber::operator*=(const ComplexNumber& number)
	{
		*this = (*this)*number;
		return *this;
	}
	ComplexNumber ComplexNumber::operator*(const double& factor) const{return ComplexNumber(real*factor,imaginary*factor);}
	ComplexNumber ComplexNumber::operator*=(const double& factor)
	{
		real *= factor;
		imaginary *= factor;
		return *this;
	}
	ComplexNumber ComplexNumber::operator/(const ComplexNumber& number) const{return (((*this)*(number.Conjugate()))/number.SquaredAmplitude());}
	ComplexNumber ComplexNumber::operator/=(const ComplexNumber& number)
	{
		*this = (*this)/number;
		return *this;
	}
	ComplexNumber ComplexNumber::operator/(const double& factor) const{return ComplexNumber(real/factor,imaginary/factor);}
	ComplexNumber ComplexNumber::operator/=(const double& factor)
	{
		real /= factor;
		imaginary /= factor;
		return *this;
	}
	ComplexNumber ComplexNumber::Conjugate() const{return ComplexNumber(real,-imaginary);}
	void ComplexNumber::Initialize()
	{
		real = 0.0;
		imaginary = 0.0;
	}
}




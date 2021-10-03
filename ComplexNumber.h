// ComplexNumber.h
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

#ifndef CASDEV_COMPLEXNUMBER_H_
#define CASDEV_COMPLEXNUMBER_H_

namespace CASDev
{
	class ComplexNumber
	{
	public:
		ComplexNumber();
		ComplexNumber(const ComplexNumber& number);
		ComplexNumber(const double& real_part,const double& imaginary_part);
		~ComplexNumber();
		ComplexNumber& operator=(const ComplexNumber& number);
		void Reset();
		void Set(const double& real_part,const double& imaginary_part);
		void AmplitudePhase(const double& amplitude,const double& phase);
		void Real(const double& value);
		void Imaginary(const double& value);
		double Real() const;
		double Imaginary() const;
		double Amplitude() const;
		double Phase() const;
		double SquaredAmplitude() const;
		bool IsReal(const double& tolerance = 1.0e-6) const;
		bool IsImaginary(const double& tolerance = 1.0e-6) const;
		ComplexNumber operator+(const ComplexNumber& number) const;
		ComplexNumber operator+=(const ComplexNumber& number);
		ComplexNumber operator-(const ComplexNumber& number) const;
		ComplexNumber operator-=(const ComplexNumber& number);
		ComplexNumber operator+(const double& number) const;
		ComplexNumber operator+=(const double& number);
		ComplexNumber operator-(const double& number) const;
		ComplexNumber operator-=(const double& number);
		ComplexNumber operator*(const ComplexNumber& number) const;
		ComplexNumber operator*=(const ComplexNumber& number);
		ComplexNumber operator*(const double& factor) const;
		ComplexNumber operator*=(const double& factor);
		ComplexNumber operator/(const ComplexNumber& number) const;
		ComplexNumber operator/=(const ComplexNumber& number);
		ComplexNumber operator/(const double& factor) const;
		ComplexNumber operator/=(const double& factor);
		ComplexNumber Conjugate() const;
		
	private:
		void Initialize();
		double real;
		double imaginary;
	};
}

#endif


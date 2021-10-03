// Signal.h 
// Copyright (c) 2014 Ahmed M. Hussein (amhussein4@gmail.com)
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

#ifndef CASDEV_SIGNAL_H_
#define CASDEV_SIGNAL_H_

#include "ComplexNumber.h"
#include "Image.h"

namespace CASDev
{
	unsigned int ReverseBits(const unsigned int& number,const unsigned int& bit_count);
	void DFT(const ComplexNumber* signal,ComplexNumber* transform,const unsigned int& size,const bool& inverse = false);
	void FFT(const ComplexNumber* signal,ComplexNumber* transform,const unsigned int& size,const bool& inverse = false);

	class Signal1D
	{
	public:
		Signal1D();
		Signal1D(const unsigned int& signal_length);
		Signal1D(const Signal1D& signal);
		~Signal1D();
		Signal1D& operator=(const Signal1D& signal);
		void Reset();
		void Allocate(const unsigned int& signal_length);
		unsigned int Length() const;
		unsigned int ActualLength() const;
		ComplexNumber operator()(const unsigned int& index) const;
		void operator()(const unsigned int& index,const ComplexNumber& value);
		void Randomize(const double& min,const double& max);
		void RandomizeReal(const double& min,const double& max);
		void RandomizeImaginary(const double& min,const double& max);
		void Increment(const unsigned int& index,const ComplexNumber& number);
		void Decrement(const unsigned int& index,const ComplexNumber& number);
		void operator*(const double& factor);
		void operator/(const double& factor);
		Signal1D Fourier() const;
		Signal1D InverseFourier() const;
		void ShiftCenter();
		void AutoCorrelation(Signal1D& autocorrelation,const bool& shift_center = true) const;
		static double Test(const unsigned int& size);

	private:
		void Initialize();
		Signal1D FourierAndInverse(const bool& inverse) const;
		unsigned int actual_length;
		unsigned int length;
		ComplexNumber* data;
		static bool use_fft;
	};

	class Signal2D
	{
	public:
		Signal2D();
		Signal2D(const unsigned int& x_size,const unsigned int& y_size);
		Signal2D(const Signal2D& signal);
		~Signal2D();
		Signal2D& operator=(const Signal2D& signal);
		void Reset();
		void Allocate(const unsigned int& x_size,const unsigned int& y_size);
		unsigned int XLength() const;
		unsigned int YLength() const;
		unsigned int ActualXLength() const;
		unsigned int ActualYLength() const;
		ComplexNumber operator()(const unsigned int& x_index,const unsigned int& y_index) const;
		void operator()(const unsigned int& x_index,const unsigned int& y_index,const ComplexNumber& value);
		void Randomize(const double& min,const double& max);
		void RandomizeReal(const double& min,const double& max);
		void RandomizeImaginary(const double& min,const double& max);
		void ReadFromPNG(const PNG& source,const char& channel);
		void WriteToPNG(PNG& target,const char& channel,const int& data_type = 0) const;
		void WriteLogAmplitudeToPNG(PNG& target) const;
		void Increment(const unsigned int& x_index,const unsigned int& y_index,const ComplexNumber& number);
		void Decrement(const unsigned int& x_index,const unsigned int& y_index,const ComplexNumber& number);
		void operator*(const double& factor);
		void operator/(const double& factor);
		void Transpose();
		Signal2D Fourier() const;
		Signal2D InverseFourier() const;
		//unsigned int PowerSpectrum(double*& spectrum) const;
		void ShiftCenter();
		void AutoCorrelation(Signal2D& autocorrelation,const bool& shift_center = true) const;
		double RadialPowerDistribution(const unsigned int& resolution,double*& distribution) const;
		double AngularPowerDistribution(const unsigned int& resolution,double*& distribution) const;
		//void SumAtProcess(const int& process_id = 0);
		static double Test(const unsigned int& x_size,const unsigned int& y_size);
		void WriteVTK(const char* filename,const double& width,const double& height) const;

	private:
		void Initialize();
		void ClearData();
		Signal2D FourierAndInverse(const bool& inverse) const;
		unsigned int actual_x_length;
		unsigned int actual_y_length;
		unsigned int x_length;
		unsigned int y_length;
		ComplexNumber** data;
		static bool use_fft;
	};

	class Signal3D
	{
	public:
		Signal3D();
		Signal3D(const unsigned int& x_size,const unsigned int& y_size,const unsigned int& z_size);
		Signal3D(const Signal3D& signal);
		~Signal3D();
		Signal3D& operator=(const Signal3D& signal);
		void Reset();
		void Allocate(const unsigned int& x_size,const unsigned int& y_size,const unsigned int& z_size);
		unsigned int XLength() const;
		unsigned int YLength() const;
		unsigned int ZLength() const;
		unsigned int ActualXLength() const;
		unsigned int ActualYLength() const;
		unsigned int ActualZLength() const;
		ComplexNumber operator()(const unsigned int& x_index,const unsigned int& y_index,const unsigned int& z_index) const;
		void operator()(const unsigned int& x_index,const unsigned int& y_index,const unsigned int& z_index,const ComplexNumber& value);
		void Randomize(const double& min,const double& max);
		void RandomizeReal(const double& min,const double& max);
		void RandomizeImaginary(const double& min,const double& max);
		//void WriteLogAmplitudeToPNG(PNG& target) const;
		void Increment(const unsigned int& x_index,const unsigned int& y_index,const unsigned int& z_index,const ComplexNumber& number);
		void Decrement(const unsigned int& x_index,const unsigned int& y_index,const unsigned int& z_index,const ComplexNumber& number);
		void operator*(const double& factor);
		void operator/(const double& factor);
		//void Transpose();
		Signal3D Fourier() const;
		Signal3D InverseFourier() const;
		////unsigned int PowerSpectrum(double*& spectrum) const;
		void ShiftCenter();
		void AutoCorrelation(Signal3D& autocorrelation,const bool& shift_center = true) const;
		//double RadialPowerDistribution(const unsigned int& resolution,double*& distribution) const;
		//double AngularPowerDistribution(const unsigned int& resolution,double*& distribution) const;
		////void SumAtProcess(const int& process_id = 0);
		static double Test(const unsigned int& x_size,const unsigned int& y_size,const unsigned int& z_size);
		//void WriteVTK(const char* filename,const double& width,const double& height) const;

	private:
		void Initialize();
		void ClearData();
		Signal3D FourierAndInverse(const bool& inverse) const;
		unsigned int actual_x_length;
		unsigned int actual_y_length;
		unsigned int actual_z_length;
		unsigned int x_length;
		unsigned int y_length;
		unsigned int z_length;
		ComplexNumber*** data;
		static bool use_fft;
	};
}

#endif


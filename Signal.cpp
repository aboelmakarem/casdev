// Signal.cpp
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

#include "CASDev.h"
#include "Signal.h"
#include "Randomizer.h"
#include "stdio.h"
#include "math.h"
#include "float.h"

namespace CASDev
{
	unsigned int ReverseBits(const unsigned int& number,const unsigned int& bit_count)
	{
		unsigned int input_mask = 1;
		unsigned int output_mask = pow(2,bit_count);
		unsigned int bit = 0;
		unsigned int result = 0;
		while(true)
		{
			// shift to the next bit to set (move left)
			output_mask >>= 1;
			if(input_mask > number)		break;
			// see if the current bit is 1
			bit = number & input_mask;
			// if it is set the corresponding but on the other end of the number to 1
			if(bit > 0)		result = result | output_mask;
			// shift to the next bit to read (move right)
			input_mask <<= 1;	
		}
		return result;
	}
	void FFT(const ComplexNumber* signal,ComplexNumber* transform,const unsigned int& size,const bool& inverse)
	{
		// this function assumes that both the signal and transform arrays have been 
		// allocated and their common size is a power of two
		// permute the input array members and set in the transform
		// the bit reversal re-orders the data so that they follow the order of the 
		// lowermost recursion level
		unsigned int depth = (unsigned int)ceil(log(size)/log(2.0));
		for(unsigned int i = 0 ; i < size ; i++)
		{
			transform[i] = signal[ReverseBits(i,depth)];
		}
		unsigned int chunk_size = 1;
		unsigned int chunk_count = size;
		// stride is the jump between the two transform entries that will be computed from the 
		// same even and odd sub-signal transformations
		unsigned int stride = 1;
		// the offset is the start of the signal chunk being processed and it is incremented by 
		// the chunk size after every chunk processing
		unsigned int offset = 0;
		double prefactor = -PI;
		if(inverse)		prefactor = PI;
		ComplexNumber factor;
		ComplexNumber term;
		for(unsigned int level = 1 ; level <= depth ; level++)
		{
			// update the chunk size of this level (multiply by 2)
			chunk_size <<= 1;
			// update the chunk count of this level (divide by 2)
			chunk_count >>= 1;
			// in every level, go over all the chunks
			offset = 0;
			for(unsigned int i = 0 ; i < chunk_count ; i++)
			{
				// go over all the entries of this chunk
				// we only need to iterate up to the stride not to the full chunk size 
				// because two transform chunk entries will be handled per iteration and the 
				// stride is always half of the chunk size
				for(unsigned int j = 0 ; j < stride ; j++)
				{
					// compute the twiddle factor, no need to multiply PI by 2 because we are 
					// dividing by the stride not by the full chunk size
					factor.AmplitudePhase(1.0,prefactor/stride*(double)j);
					// multiply the twiddle factor by the term from the lower level odd 
					// transform
					term = factor*transform[offset + j + stride];
					// set the term after the stride
					transform[offset + j + stride] = transform[offset + j] - term;
					// set the term before the stride
					transform[offset + j] = transform[offset + j] + term;
				}
				// move to the next chunk
				offset += chunk_size;
			}
			// double the size of the stride, the stride has to always be equal to 
			// half the signal size at this level
			stride = chunk_size;
		}
		if(inverse)
		{
			double scale = (double)size;
			for(unsigned int i = 0 ; i < size ; i++)
			{
				transform[i] /= scale;
			}
		}
	}
	void DFT(const ComplexNumber* signal,ComplexNumber* transform,const unsigned int& size,const bool& inverse)
	{
		// this function assumes that both the function and transform arrays have been 
		// allocated
		double prefactor = -2.0*PI/(double)size;
		if(inverse)			prefactor *= -1.0;
		ComplexNumber factor;
		for(unsigned int i = 0 ; i < size ; i++)
		{
			transform[i].Reset();
			for(unsigned int j = 0 ; j < size ; j++)
			{
				factor.AmplitudePhase(1.0,prefactor*(double)i*(double)j);
				transform[i] += factor*signal[j];
			}
		}
		if(inverse)
		{
			double scale = (double)size;
			for(unsigned int i = 0 ; i < size ; i++)
			{
				transform[i] /= scale;
			}
		}
	}

	bool Signal1D::use_fft = true;
	Signal1D::Signal1D(){Initialize();}
	Signal1D::Signal1D(const unsigned int& signal_length)
	{
		Initialize();
		Allocate(signal_length);
	}
	Signal1D::Signal1D(const Signal1D& signal)
	{
		Initialize();
		*this = signal;
	}
	Signal1D::~Signal1D(){Reset();}
	Signal1D& Signal1D::operator=(const Signal1D& signal)
	{
		Allocate(signal.length);
		for(unsigned int i = 0 ; i < length ; i++)
		{
			data[i] = signal.data[i];
		}
		actual_length = signal.actual_length;
		return *this;
	}
	void Signal1D::Reset()
	{
		if(data != 0)		delete [] data;
		Initialize();
	}
	void Signal1D::Allocate(const unsigned int& signal_length)
	{
		if(length == signal_length)				return;
		// allocate a length that is a power of two
		actual_length = signal_length;

		double real_exponent = log(actual_length)/log(2.0);
		double tolerance = 1.0e-6;
		if(fabs(real_exponent - floor(real_exponent)) < tolerance)		length = actual_length;
		else
		{
			int exponent = (int)floor(real_exponent) + 1;
			length = (unsigned int)floor(pow(2.0,exponent) + 0.5);
		}
		if(data != 0)			delete [] data;
		data = new ComplexNumber[length];
		for(unsigned int i = 0 ; i < length ; i++)
		{
			data[i].Set(0.0,0.0);
		}
	}
	unsigned int Signal1D::Length() const{return length;}
	unsigned int Signal1D::ActualLength() const{return actual_length;}
	ComplexNumber Signal1D::operator()(const unsigned int& index) const{return data[index];}
	void Signal1D::operator()(const unsigned int& index,const ComplexNumber& value){data[index] = value;}
	void Signal1D::Randomize(const double& min,const double& max)
	{
		for(unsigned int i = 0 ; i < actual_length ; i++)
		{
			data[i].Set(Randomizer::Uniform(min,max),Randomizer::Uniform(min,max));
		}
	}
	void Signal1D::RandomizeReal(const double& min,const double& max)
	{
		for(unsigned int i = 0 ; i < actual_length ; i++)
		{
			data[i].Real(Randomizer::Uniform(min,max));
		}
	}
	void Signal1D::RandomizeImaginary(const double& min,const double& max)
	{
		for(unsigned int i = 0 ; i < actual_length ; i++)
		{
			data[i].Imaginary(Randomizer::Uniform(min,max));
		}
	}
	void Signal1D::Increment(const unsigned int& index,const ComplexNumber& number){data[index] += number;}
	void Signal1D::Decrement(const unsigned int& index,const ComplexNumber& number){data[index] -= number;}
	void Signal1D::operator*(const double& factor)
	{
		for(unsigned int i = 0 ; i < length ; i++)
		{
			data[i] *= factor;
		}
	}
	void Signal1D::operator/(const double& factor)
	{
		for(unsigned int i = 0 ; i < length ; i++)
		{
			data[i] /= factor;
		}
	}
	Signal1D Signal1D::Fourier() const{return FourierAndInverse(false);}
	Signal1D Signal1D::InverseFourier() const{return FourierAndInverse(true);}
	void Signal1D::ShiftCenter()
	{
		ComplexNumber* shifted_data = new ComplexNumber[length];
		unsigned int shift = length/2;
		for(unsigned int i = 0 ; i < length ; i++)
		{
			shifted_data[(i + shift)%length] = data[i];
		}
		delete [] data;
		data = shifted_data;
	}
	void Signal1D::AutoCorrelation(Signal1D& autocorrelation,const bool& shift_center) const
	{
		// the Fourier transform of a signal's autocorrelation 
		// is the product of the conjugate of its Fourier transform 
		// by its Fourier transform
		Signal1D transform = Fourier();
		Signal1D product(length);
		for(unsigned int i = 0 ; i < length ; i++)
		{
			product(i,ComplexNumber(transform(i).SquaredAmplitude(),0.0));
		}
		// inverse tranform the product to get the autocorrelation
		autocorrelation = product.InverseFourier();
		// shift its center if required
		if(shift_center)		autocorrelation.ShiftCenter();
	}
	double Signal1D::Test(const unsigned int& size)
	{
		Signal1D signal(size);
		signal.Randomize(-100.0,100.0);
		//printf("signal:\n");
		//for(unsigned int i = 0 ; i < size ; i++)
		//{
		//	printf("\t%f + %f i\n",signal(i).Real(),signal(i).Imaginary());
		//}
		use_fft = false;
		//Signal1D direct_transform = signal.Fourier();
		use_fft = true;
		Signal1D fft_transform = signal.Fourier();
		//printf("transform:\n");
		//for(unsigned int i = 0 ; i < size ; i++)
		//{
		//	printf("\t%f + %f i : %f + %f i\n",direct_transform(i).Real(),direct_transform(i).Imaginary(),fft_transform(i).Real(),fft_transform(i).Imaginary());
		//}
		use_fft = false;
		//Signal1D direct_direct_recovery = direct_transform.InverseFourier();
		//Signal1D direct_fft_recovery = fft_transform.InverseFourier();
		use_fft = true;
		//Signal1D fft_direct_recovery = direct_transform.InverseFourier();
		Signal1D fft_fft_recovery = fft_transform.InverseFourier();
		//printf("recovered signal:\n");
		//for(unsigned int i = 0 ; i < size ; i++)
		//{
		//	printf("DD: \t%f + %f i : %f + %f i\n",signal(i).Real(),signal(i).Imaginary(),direct_direct_recovery(i).Real(),direct_direct_recovery(i).Imaginary());
		//	printf("DF: \t%f + %f i : %f + %f i\n",signal(i).Real(),signal(i).Imaginary(),fft_direct_recovery(i).Real(),fft_direct_recovery(i).Imaginary());
		//	printf("FD: \t%f + %f i : %f + %f i\n",signal(i).Real(),signal(i).Imaginary(),direct_fft_recovery(i).Real(),direct_fft_recovery(i).Imaginary());
		//	printf("FF: \t%f + %f i : %f + %f i\n",signal(i).Real(),signal(i).Imaginary(),fft_fft_recovery(i).Real(),fft_fft_recovery(i).Imaginary());
		//}
		ComplexNumber error;
		//double direct_direct_rms_error = 0.0;
		//double fft_direct_rms_error = 0.0;
		//double direct_fft_rms_error = 0.0;
		double fft_fft_rms_error = 0.0;
		for(unsigned int i = 0 ; i < size ; i++)
		{
			//error = signal(i) - direct_direct_recovery(i);
			//direct_direct_rms_error += error.SquaredAmplitude();
			//error = signal(i) - fft_direct_recovery(i);
			//fft_direct_rms_error += error.SquaredAmplitude();
			//error = signal(i) - direct_fft_recovery(i);
			//direct_fft_rms_error += error.SquaredAmplitude();
			error = signal(i) - fft_fft_recovery(i);
			fft_fft_rms_error += error.SquaredAmplitude();
		}
		//direct_direct_rms_error = sqrt(direct_direct_rms_error)/(double)size;
		//fft_direct_rms_error = sqrt(fft_direct_rms_error)/(double)size;
		//direct_fft_rms_error = sqrt(direct_fft_rms_error)/(double)size;
		fft_fft_rms_error = sqrt(fft_fft_rms_error/(double)size);
		//printf("direct-direct recovery rms error : %e\n",direct_direct_rms_error);
		//printf("fft-direct recovery rms error : %e\n",fft_direct_rms_error);
		//printf("direct-fft recovery rms error : %e\n",direct_fft_rms_error);
		printf("fft-fft recovery rms error : %e\n",fft_fft_rms_error);
		return fft_fft_rms_error;
	}
	void Signal1D::Initialize()
	{
		actual_length = 0;
		length = 0;
		data = 0;
	}
	Signal1D Signal1D::FourierAndInverse(const bool& inverse) const
	{
		Signal1D transform(length);
		if(use_fft)			FFT(data,transform.data,length,inverse);
		else				DFT(data,transform.data,length,inverse);
		transform.actual_length = actual_length;
		return transform;
	}

	bool Signal2D::use_fft = true;
	Signal2D::Signal2D(){Initialize();}
	Signal2D::Signal2D(const unsigned int& x_size,const unsigned int& y_size)
	{
		Initialize();
		Allocate(x_size,y_size);
	}
	Signal2D::Signal2D(const Signal2D& signal)
	{
		Initialize();
		*this = signal;
	}
	Signal2D::~Signal2D(){Reset();}
	Signal2D& Signal2D::operator=(const Signal2D& signal)
	{
		Allocate(signal.x_length,signal.y_length);
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				data[i][j] = signal.data[i][j];
			}
		}
		actual_x_length = signal.actual_x_length;
		actual_y_length = signal.actual_y_length;
		return *this;
	}
	void Signal2D::Reset()
	{
		ClearData();
		Initialize();
	}
	void Signal2D::Allocate(const unsigned int& x_size,const unsigned int& y_size)
	{
		if((x_length == x_size) && (y_length == y_size))						return;
		// if the fast fourier transform will be used, then allocate a length 
		// that is a power of two
		// allocate a length that is a power of two for both dimensions
		actual_x_length = x_size;
		actual_y_length = y_size;
		double tolerance = 1.0e-6;
		int exponent = 0;
		// x dimension length
		double real_exponent = log(actual_x_length)/log(2.0);
		if(fabs(real_exponent - floor(real_exponent)) < tolerance)		x_length = actual_x_length;
		else
		{
			exponent = (int)floor(real_exponent) + 1;
			x_length = (unsigned int)floor(pow(2.0,exponent) + 0.5);
		}
		// y dimension length
		real_exponent = log(actual_y_length)/log(2.0);
		if(fabs(real_exponent - floor(real_exponent)) < tolerance)		y_length = actual_y_length;
		else
		{
			exponent = (int)floor(real_exponent) + 1;
			y_length = (unsigned int)floor(pow(2.0,exponent) + 0.5);
		}
		ClearData();
		data = new ComplexNumber*[x_length];
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			data[i] = new ComplexNumber[y_length];
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				data[i][j].Set(0.0,0.0);
			}
		}
	}
	unsigned int Signal2D::XLength() const{return x_length;}
	unsigned int Signal2D::YLength() const{return y_length;}
	unsigned int Signal2D::ActualXLength() const{return actual_x_length;}
	unsigned int Signal2D::ActualYLength() const{return actual_y_length;}
	ComplexNumber Signal2D::operator()(const unsigned int& x_index,const unsigned int& y_index) const{return data[x_index][y_index];}
	void Signal2D::operator()(const unsigned int& x_index,const unsigned int& y_index,const ComplexNumber& value){data[x_index][y_index] = value;}
	void Signal2D::Randomize(const double& min,const double& max)
	{
		for(unsigned int i = 0 ; i < actual_x_length ; i++)
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				data[i][j].Set(Randomizer::Uniform(min,max),Randomizer::Uniform(min,max));
			}
		}
	}
	void Signal2D::RandomizeReal(const double& min,const double& max)
	{
		for(unsigned int i = 0 ; i < actual_x_length ; i++)
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				data[i][j].Real(Randomizer::Uniform(min,max));
			}
		}
	}
	void Signal2D::RandomizeImaginary(const double& min,const double& max)
	{
		for(unsigned int i = 0 ; i < actual_x_length ; i++)
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				data[i][j].Imaginary(Randomizer::Uniform(min,max));
			}
		}
	}
	void Signal2D::ReadFromPNG(const PNG& source,const char& channel)
	{
		Reset();
		unsigned int width = source.Width();
		unsigned int height = source.Height();
		Allocate(width,height);
		double factor = 1.0/255.0;
		if(channel == 'R')
		{
			for(unsigned int j = 0 ; j < height ; j++)
			{
				for(unsigned int i = 0 ; i < width ; i++)
				{
					data[i][j] = ComplexNumber(factor*((double)source.R(j,i)),0.0);
				}
			}
		}
		if(channel == 'G')
		{
			for(unsigned int j = 0 ; j < height ; j++)
			{
				for(unsigned int i = 0 ; i < width ; i++)
				{
					data[i][j] = ComplexNumber(factor*((double)source.G(j,i)),0.0);
				}
			}
		}
		if(channel == 'B')
		{
			for(unsigned int j = 0 ; j < height ; j++)
			{
				for(unsigned int i = 0 ; i < width ; i++)
				{
					data[i][j] = ComplexNumber(factor*((double)source.B(j,i)),0.0);
				}
			}
		}
		if(channel == 'A')
		{
			for(unsigned int j = 0 ; j < height ; j++)
			{
				for(unsigned int i = 0 ; i < width ; i++)
				{
					data[i][j] = ComplexNumber(factor*((double)source.A(j,i)),0.0);
				}
			}
		}
	}
	void Signal2D::WriteToPNG(PNG& target,const char& channel,const int& data_type) const
	{
		// data type is the data to write to the PNG channel
		// data_type = 0 : write signal real data
		// data_type = 1 : write signal imaginary data
		// data_type = 2 : write signal amplitude data
		// data_type = 3 : write signal phase data
		int write_type = 0;
		if(data_type == 1)		write_type = 1;
		else if(data_type == 2)		write_type = 2;
		else if(data_type == 3)		write_type = 3;
		double min = DBL_MAX;
		double max = -DBL_MAX;
		double value = 0.0;
		// get the min and max values
		for(unsigned int j = 0 ; j < actual_y_length ; j++)
		{
			for(unsigned int i = 0 ; i < actual_x_length ; i++)
			{
				if(write_type == 0)			value = data[i][j].Real();
				else if(write_type == 1)	value = data[i][j].Imaginary();
				else if(write_type == 2)	value = data[i][j].Amplitude();
				else if(write_type == 3)	value = data[i][j].Phase();
				else if(write_type == 4)	value = log(data[i][j].Amplitude());
				if(value < min)			min = value;
				else if(value > max)	max = value;
			}
		}
		double factor = 255.0/(max - min);
		if(channel == 'R')
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				for(unsigned int i = 0 ; i < actual_x_length ; i++)
				{
					if(write_type == 0)				value = factor*(data[i][j].Real() - min);
					else if(write_type == 1)		value = factor*(data[i][j].Imaginary() - min);
					else if(write_type == 2)		value = factor*(data[i][j].Amplitude() - min);
					else if(write_type == 3)		value = factor*(data[i][j].Phase() - min);
					else if(write_type == 4)		value = factor*(log(data[i][j].Amplitude()) - min);
					target.R(j,i,(unsigned char)floor(value + 0.5));
				}
			}
		}
		else if(channel == 'G')
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				for(unsigned int i = 0 ; i < actual_x_length ; i++)
				{
					if(write_type == 0)				value = factor*(data[i][j].Real() - min);
					else if(write_type == 1)		value = factor*(data[i][j].Imaginary() - min);
					else if(write_type == 2)		value = factor*(data[i][j].Amplitude() - min);
					else if(write_type == 3)		value = factor*(data[i][j].Phase() - min);
					else if(write_type == 4)		value = factor*(log(data[i][j].Amplitude()) - min);
					target.G(j,i,(unsigned char)floor(value + 0.5));
				}
			}
		}
		else if(channel == 'B')
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				for(unsigned int i = 0 ; i < actual_x_length ; i++)
				{
					if(write_type == 0)				value = factor*(data[i][j].Real() - min);
					else if(write_type == 1)		value = factor*(data[i][j].Imaginary() - min);
					else if(write_type == 2)		value = factor*(data[i][j].Amplitude() - min);
					else if(write_type == 3)		value = factor*(data[i][j].Phase() - min);
					else if(write_type == 4)		value = factor*(log(data[i][j].Amplitude()) - min);
					target.B(j,i,(unsigned char)floor(value + 0.5));
				}
			}
		}
		else if(channel == 'A')
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				for(unsigned int i = 0 ; i < actual_x_length ; i++)
				{
					if(write_type == 0)				value = factor*(data[i][j].Real() - min);
					else if(write_type == 1)		value = factor*(data[i][j].Imaginary() - min);
					else if(write_type == 2)		value = factor*(data[i][j].Amplitude() - min);
					else if(write_type == 3)		value = factor*(data[i][j].Phase() - min);
					else if(write_type == 4)		value = factor*(log(data[i][j].Amplitude()) - min);
					target.A(j,i,(unsigned char)floor(value + 0.5));
				}
			}
		}
	}
	void Signal2D::WriteLogAmplitudeToPNG(PNG& target) const
	{
		// the PNG has to be preallocated for the correct size at least 3 channels
		double min = DBL_MAX;
		double max = -DBL_MAX;
		double value = 0.0;
		// get the min and max values
		double tolerance = 1.0e-12;
		for(unsigned int j = 0 ; j < y_length ; j++)
		{
			for(unsigned int i = 0 ; i < x_length ; i++)
			{
				value = data[i][j].Amplitude();
				if(value < tolerance)		value = tolerance;
				value = log(value);
				if(value < min)			min = value;
				else if(value > max)	max = value;
			}
		}
		double factor = 512.0/(max - min);
		for(unsigned int j = 0 ; j < y_length ; j++)
		{
			for(unsigned int i = 0 ; i < x_length ; i++)
			{
				value = data[i][j].Amplitude();
				if(value < tolerance)		value = tolerance;
				value = factor*(log(value) - min);
				if(value < 256.0)
				{
					target.R(j,i,0);
					target.G(j,i,(unsigned char)floor(255.0/256.0*value + 0.5));
					target.B(j,i,(unsigned char)floor(255.0*(1.0 - value/256.0) + 0.5));
				}
				else
				{
					target.R(j,i,(unsigned char)floor(255.0/256.0*value + 0.5));
					target.G(j,i,(unsigned char)floor(255.0*(1.0 - value/256.0) + 0.5));
					target.B(j,i,0);
				}
			}
		}
	}
	void Signal2D::Increment(const unsigned int& x_index,const unsigned int& y_index,const ComplexNumber& number){data[x_index][y_index] += number;}
	void Signal2D::Decrement(const unsigned int& x_index,const unsigned int& y_index,const ComplexNumber& number){data[x_index][y_index] -= number;}
	void Signal2D::operator*(const double& factor)
	{
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				data[i][j] *= factor;
			}
		}
	}
	void Signal2D::operator/(const double& factor)
	{
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				data[i][j] /= factor;
			}
		}
	}
	void Signal2D::Transpose()
	{
		ComplexNumber** transposed_data = new ComplexNumber*[y_length];
		for(unsigned int i = 0 ; i < y_length ; i++)
		{
			transposed_data[i] = new ComplexNumber[x_length];
		}
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				transposed_data[j][i] = data[i][j];
			}
			delete [] data[i];
		}
		delete [] data;
		data = transposed_data;
		unsigned int temp = x_length;
		x_length = y_length;
		y_length = temp;
		temp = actual_x_length;
		actual_x_length = actual_y_length;
		actual_y_length = temp;
	}
	Signal2D Signal2D::Fourier() const{return FourierAndInverse(false);}
	Signal2D Signal2D::InverseFourier() const{return FourierAndInverse(true);}
	/*unsigned int Signal2D::PowerSpectrum(double*& spectrum) const
	{
		// get the two-dimensional power array
		double** power_array = new double*[x_length];
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			power_array[i] = new double[y_length];
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				power_array[i][j] = data[i][j].SquaredAmplitude();
			}
		}
		// arrange into a one-dimensional power spectrum array
		unsigned int spectrum_size = (unsigned int)ceil(sqrt(x_length*x_length + y_length*y_length));
		spectrum = new double[spectrum_size];
		unsigned int* counts = new unsigned int[spectrum_size];
		for(unsigned int i = 0 ; i < spectrum_size ; i++)
		{
			spectrum[i] = 0.0;
			counts[i] = 0;
		}
		unsigned int index = 0;
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{	
				index = (unsigned int)floor(sqrt(i*i + j*j) + 0.5);
				spectrum[index] += power_array[i][j];
				counts[index]++;
			}
			delete [] power_array[i];
		}
		delete [] power_array;
		// get the average value for each index in the spectrum array
		for(unsigned int i = 0 ; i < spectrum_size ; i++)
		{
			if(counts[i] == 0)			continue;
			spectrum[i] = spectrum[i]/(double)counts[i];
		}
		delete [] counts;
		return spectrum_size;
	}*/
	void Signal2D::ShiftCenter()
	{
		ComplexNumber** shifted_data = new ComplexNumber*[x_length];
		unsigned int x_shift = x_length/2;
		unsigned int y_shift = y_length/2;
		unsigned int i_index = 0;
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			i_index = (i + x_shift)%x_length;
			shifted_data[i_index] = new ComplexNumber[y_length];
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				shifted_data[i_index][(j + y_shift)%y_length] = data[i][j];
			}
		}
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			delete [] data[i];
		}
		delete [] data;
		data = shifted_data;
	}
	void Signal2D::AutoCorrelation(Signal2D& autocorrelation,const bool& shift_center) const
	{
		// the Fourier transform of a signal's autocorrelation 
		// is the product of the conjugate of its Fourier transform 
		// by its Fourier transform
		Signal2D transform = Fourier();
		Signal2D product(x_length,y_length);
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				product(i,j,ComplexNumber(transform(i,j).SquaredAmplitude(),0.0));
			}
		}
		// inverse tranform the product to get the autocorrelation
		autocorrelation = product.InverseFourier();
		// shift its center if required
		if(shift_center)		autocorrelation.ShiftCenter();
	}
	double Signal2D::RadialPowerDistribution(const unsigned int& resolution,double*& distribution) const
	{
		// this function assumes that the signal is center shifted
		distribution = new double[resolution];
		unsigned int* counts = new unsigned int[resolution];
		for(unsigned int i = 0 ; i < resolution ; i++)
		{
			distribution[i] = 0.0;
			counts[i] = 0;
		}
		// 0.5 factor because the signal is center shifted and max i/j are at half x/y lengths
		double x = 0.5*(double)x_length;
		double y = 0.5*(double)y_length;
		double delta = sqrt(x*x + y*y)/(double)resolution;
		unsigned int index = 0;
		double x_shift = 0.5*(double)x_length;
		double y_shift = 0.5*(double)y_length;
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{	
				x = double(i) - x_shift;
				y = double(j) - y_shift;
				index = (unsigned int)floor(sqrt(x*x + y*y)/delta);
				distribution[index] += data[i][j].SquaredAmplitude();
				counts[index]++;
			}
		}
		// get the average value for each bin in the distribution array
		for(unsigned int i = 0 ; i < resolution ; i++)
		{
			if(counts[i] == 0)			continue;
			distribution[i] = distribution[i]/(double)counts[i];
		}
		delete [] counts;
		return delta;
	}
	double Signal2D::AngularPowerDistribution(const unsigned int& resolution,double*& distribution) const
	{
		// this function assumes that the signal is center shifted
		distribution = new double[resolution];
		unsigned int* counts = new unsigned int[resolution];
		for(unsigned int i = 0 ; i < resolution ; i++)
		{
			distribution[i] = 0.0;
			counts[i] = 0;
		}
		double delta = 2.0*PI/(double)resolution;
		unsigned int index = 0;
		double x = 0.0;
		double y = 0.0;
		double x_shift = 0.5*(double)x_length;
		double y_shift = 0.5*(double)y_length;
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{	
				// skip the origin
				if((i == x_length/2) && (j == y_length/2))		continue;
				x = double(i) - x_shift;
				y = double(j) - y_shift;
				index = (unsigned int)floor((atan2(y,x) + PI)/delta);
				distribution[index] += data[i][j].SquaredAmplitude();
				counts[index]++;
			}
		}
		// get the average value for each bin in the distribution array
		for(unsigned int i = 0 ; i < resolution ; i++)
		{
			if(counts[i] == 0)			continue;
			distribution[i] = distribution[i]/(double)counts[i];
		}
		delete [] counts;
		return delta;
	}
	/*void Complex2DArray::SumAtProcess(const int& process_id)
	{
		unsigned int entry_count = 2*x_length*y_length;
		double* entries = new double[entry_count];
		double* sums = 0;
		if(EZ::Parallel::EZMPI::ProcessID() == process_id)		sums = new double[entry_count];
		unsigned int index = 0;
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				entries[index++] = data[i][j].Real();
				entries[index++] = data[i][j].Imaginary();
			}
		}
		EZ::Parallel::EZMPI::Reduce(entries,sums,entry_count,MPI_DOUBLE,MPI_SUM,process_id);
		delete [] entries;
		EZ::Parallel::EZMPI::Barrier();
		if(EZ::Parallel::EZMPI::ProcessID() == process_id)
		{
			index = 0;
			for(unsigned int i = 0 ; i < x_length ; i++)
			{
				for(unsigned int j = 0 ; j < y_length ; j++)
				{
					data[i][j].Real(sums[index++]);
					data[i][j].Imaginary(sums[index++]);
				}
			}
			delete [] sums;
		}
		EZ::Parallel::EZMPI::Barrier();
	}*/
	double Signal2D::Test(const unsigned int& x_size,const unsigned int& y_size)
	{
		Signal2D signal(x_size,y_size);
		signal.Randomize(-100.0,100.0);
		//printf("signal:\n");
		//for(unsigned int i = 0 ; i < x_size ; i++)
		//{
			//for(unsigned int j = 0 ; j < y_size ; j++)
			//{
				//printf("\t%f + %f i\n",signal(i,j).Real(),signal(i,j).Imaginary());
			//}
		//}
		//use_fft = false;
		//Signal2D direct_transform = signal.Fourier();
		//use_fft = true;
		Signal2D fft_transform = signal.Fourier();
		//printf("transform:\n");
		//for(unsigned int i = 0 ; i < x_size ; i++)
		//{
			//for(unsigned int j = 0 ; j < y_size ; j++)
			//{
				//printf("\t%f + %f i\n",fft_transform(i).Real(),fft_transform(i).Imaginary());
			//}
		//}
		//use_fft = false;
		//Signal2D direct_direct_recovery = direct_transform.InverseFourier();
		//Signal2D direct_fft_recovery = fft_transform.InverseFourier();
		//use_fft = true;
		//Signal2D fft_direct_recovery = direct_transform.InverseFourier();
		Signal2D fft_fft_recovery = fft_transform.InverseFourier();
		//printf("recovered signal:\n");
		//for(unsigned int i = 0 ; i < x_size ; i++)
		//{
			//for(unsigned int j = 0 ; j < y_size ; j++)
			//{
				//printf("\t%f + %f i : %f + %f i\n",signal(i).Real(),signal(i).Imaginary(),recovery(i).Real(),recovery(i).Imaginary());
			//}
		//}
		ComplexNumber error;
		//double direct_direct_rms_error = 0.0;
		//double fft_direct_rms_error = 0.0;
		//double direct_fft_rms_error = 0.0;
		double fft_fft_rms_error = 0.0;
		for(unsigned int i = 0 ; i < x_size ; i++)
		{
			for(unsigned int j = 0 ; j < y_size ; j++)
			{
				//error = signal(i) - direct_direct_recovery(i);
				//direct_direct_rms_error += error.SquaredAmplitude();
				//error = signal(i) - fft_direct_recovery(i);
				//fft_direct_rms_error += error.SquaredAmplitude();
				//error = signal(i) - direct_fft_recovery(i);
				//direct_fft_rms_error += error.SquaredAmplitude();
				error = signal(i,j) - fft_fft_recovery(i,j);
				fft_fft_rms_error += error.SquaredAmplitude();
			}
		}
		//direct_direct_rms_error = sqrt(direct_direct_rms_error)/(double)size;
		//fft_direct_rms_error = sqrt(fft_direct_rms_error)/(double)size;
		//direct_fft_rms_error = sqrt(direct_fft_rms_error)/(double)size;
		fft_fft_rms_error = sqrt(fft_fft_rms_error/(double)x_size/(double)y_size);
		//printf("direct-direct recovery rms error : %e\n",direct_direct_rms_error);
		//printf("fft-direct recovery rms error : %e\n",fft_direct_rms_error);
		//printf("direct-fft recovery rms error : %e\n",direct_fft_rms_error);
		printf("fft-fft recovery rms error : %e\n",fft_fft_rms_error);
		return fft_fft_rms_error;
	}
	/*void Signal2D::WriteVTK(const char* filename,const double& width,const double& height) const
	{
		double dx = width/(double)(actual_x_length - 1);
		double dy = height/(double)(actual_y_length - 1);
		unsigned int point_count = actual_x_length*actual_y_length;
		FILE* file = fopen(filename,"w");
		fprintf(file,"# vtk DataFile Version 3.0\n");
		fprintf(file,"Two-dimensional signal\n");
		fprintf(file,"ASCII\n");
		fprintf(file,"DATASET UNSTRUCTURED_GRID\n");
		fprintf(file,"POINTS %u double\n",point_count);
		for(unsigned int i = 0 ; i < actual_x_length ; i++)
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				fprintf(file,"%15.10f\t%15.10f\t%15.10f\n",i*dx,j*dy,sqrt(data[i][j].SquaredAmplitude()));
			}
		}

		// write the cell connectivity and their types
		unsigned int cell_count = 2*(actual_x_length - 1)*(actual_y_length - 1);
		fprintf(file,"CELLS %u %u\n",cell_count,4*cell_count);
		for(unsigned int i = 0 ; i < (actual_x_length - 1) ; i++)
		{
			for(unsigned int j = 0 ; j < (actual_y_length - 1) ; j++)
			{
				fprintf(file,"3\t%u\t%u\t%u\n",i*actual_y_length + j,(i + 1)*actual_y_length + j,i*actual_y_length + j + 1);
				fprintf(file,"3\t%u\t%u\t%u\n",(i + 1)*actual_y_length + j,(i + 1)*actual_y_length + j + 1,i*actual_y_length + j + 1);
			}
		}
		fprintf(file,"CELL_TYPES %d\n",cell_count);
		for(unsigned int i = 0 ; i < cell_count ; i++)
		{
			fprintf(file,"5\n");
		}
		// write the cell data
		fprintf(file,"POINT_DATA %d\n",point_count);
		fprintf(file,"SCALARS real double\n");
		fprintf(file,"LOOKUP_TABLE default\n");
		for(unsigned int i = 0 ; i < actual_x_length ; i++)
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				fprintf(file,"%e\n",data[i][j].Real());
			}
		}
		fprintf(file,"SCALARS imaginary double\n");
		fprintf(file,"LOOKUP_TABLE default\n");
		for(unsigned int i = 0 ; i < actual_x_length ; i++)
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				fprintf(file,"%e\n",data[i][j].Imaginary());
			}
		}
		fprintf(file,"SCALARS log_amplitude double\n");
		fprintf(file,"LOOKUP_TABLE default\n");
		for(unsigned int i = 0 ; i < actual_x_length ; i++)
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				fprintf(file,"%e\n",log(sqrt(data[i][j].SquaredAmplitude())));
			}
		}
		fclose(file);
	}*/
	void Signal2D::Initialize()
	{
		actual_x_length = 0;
		actual_y_length = 0;
		x_length = 0;
		y_length = 0;
		data = 0;
	}
	void Signal2D::ClearData()
	{
		if(data == 0)		return;
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			if(data[i] != 0)		delete [] data[i];
		}
		delete [] data;
		data = 0;
		x_length = 0;
		y_length = 0;
		actual_x_length = 0;
		actual_y_length = 0;
	}
	Signal2D Signal2D::FourierAndInverse(const bool& inverse) const
	{
		Signal2D column_transform(x_length,y_length);
		// get the one dimensional FFT transforms of signal columns
		if(use_fft)
		{
			for(unsigned int i = 0 ; i < x_length ; i++)
			{
				FFT(data[i],column_transform.data[i],y_length,inverse);
			}
		}
		else
		{
			for(unsigned int i = 0 ; i < x_length ; i++)
			{
				DFT(data[i],column_transform.data[i],y_length,inverse);
			}
		}
		// transpose the signal so that we can pass the rows to FFT
		column_transform.Transpose();
		// transform the transposed transform columns (original transform rows)
		Signal2D transform(y_length,x_length);
		if(use_fft)
		{
			for(unsigned int i = 0 ; i < y_length ; i++)
			{
				FFT(column_transform.data[i],transform.data[i],x_length,inverse);
			}
		}
		else
		{
			for(unsigned int i = 0 ; i < x_length ; i++)
			{
				DFT(column_transform.data[i],transform.data[i],x_length,inverse);
			}
		}
		// transpose the transform back
		transform.Transpose();
		transform.actual_x_length = actual_x_length;
		transform.actual_y_length = actual_y_length;
		return transform;
	}

	bool Signal3D::use_fft = true;
	Signal3D::Signal3D(){Initialize();}
	Signal3D::Signal3D(const unsigned int& x_size,const unsigned int& y_size,const unsigned int& z_size)
	{
		Initialize();
		Allocate(x_size,y_size,z_size);
	}
	Signal3D::Signal3D(const Signal3D& signal)
	{
		Initialize();
		*this = signal;
	}
	Signal3D::~Signal3D(){Reset();}
	Signal3D& Signal3D::operator=(const Signal3D& signal)
	{
		Allocate(signal.x_length,signal.y_length,signal.z_length);
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				for(unsigned int k = 0 ; k < z_length ; k++)
				{
					data[i][j][k] = signal.data[i][j][k];
				}
			}
		}
		actual_x_length = signal.actual_x_length;
		actual_y_length = signal.actual_y_length;
		actual_z_length = signal.actual_z_length;
		return *this;
	}
	void Signal3D::Reset()
	{
		ClearData();
		Initialize();
	}
	void Signal3D::Allocate(const unsigned int& x_size,const unsigned int& y_size,const unsigned int& z_size)
	{
		if((x_length == x_size) && (y_length == y_size) && (z_length == z_size))			return;
		// if the fast fourier transform will be used, then allocate a length 
		// that is a power of two
		// allocate a length that is a power of two for all dimensions
		actual_x_length = x_size;
		actual_y_length = y_size;
		actual_z_length = z_size;
		double tolerance = 1.0e-6;
		int exponent = 0;
		// x dimension length
		double real_exponent = log(actual_x_length)/log(2.0);
		if(fabs(real_exponent - floor(real_exponent)) < tolerance)		x_length = actual_x_length;
		else
		{
			exponent = (int)floor(real_exponent) + 1;
			x_length = (unsigned int)floor(pow(2.0,exponent) + 0.5);
		}
		// y dimension length
		real_exponent = log(actual_y_length)/log(2.0);
		if(fabs(real_exponent - floor(real_exponent)) < tolerance)		y_length = actual_y_length;
		else
		{
			exponent = (int)floor(real_exponent) + 1;
			y_length = (unsigned int)floor(pow(2.0,exponent) + 0.5);
		}
		// z dimension length
		real_exponent = log(actual_z_length)/log(2.0);
		if(fabs(real_exponent - floor(real_exponent)) < tolerance)		z_length = actual_z_length;
		else
		{
			exponent = (int)floor(real_exponent) + 1;
			z_length = (unsigned int)floor(pow(2.0,exponent) + 0.5);
		}
		ClearData();
		data = new ComplexNumber**[x_length];
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			data[i] = new ComplexNumber*[y_length];
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				data[i][j] = new ComplexNumber[z_length];
				for(unsigned int k = 0 ; k < z_length ; k++)
				{
					data[i][j][k].Set(0.0,0.0);
				}
			}
		}
	}
	unsigned int Signal3D::XLength() const{return x_length;}
	unsigned int Signal3D::YLength() const{return y_length;}
	unsigned int Signal3D::ZLength() const{return z_length;}
	unsigned int Signal3D::ActualXLength() const{return actual_x_length;}
	unsigned int Signal3D::ActualYLength() const{return actual_y_length;}
	unsigned int Signal3D::ActualZLength() const{return actual_z_length;}
	ComplexNumber Signal3D::operator()(const unsigned int& x_index,const unsigned int& y_index,const unsigned int& z_index) const{return data[x_index][y_index][z_index];}
	void Signal3D::operator()(const unsigned int& x_index,const unsigned int& y_index,const unsigned int& z_index,const ComplexNumber& value){data[x_index][y_index][z_index] = value;}
	void Signal3D::Randomize(const double& min,const double& max)
	{
		for(unsigned int i = 0 ; i < actual_x_length ; i++)
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				for(unsigned int k = 0 ; k < actual_z_length ; k++)
				{
					data[i][j][k].Set(Randomizer::Uniform(min,max),Randomizer::Uniform(min,max));
				}
			}
		}
	}
	void Signal3D::RandomizeReal(const double& min,const double& max)
	{
		for(unsigned int i = 0 ; i < actual_x_length ; i++)
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				for(unsigned int k = 0 ; k < actual_z_length ; k++)
				{
					data[i][j][k].Real(Randomizer::Uniform(min,max));
				}
			}
		}
	}
	void Signal3D::RandomizeImaginary(const double& min,const double& max)
	{
		for(unsigned int i = 0 ; i < actual_x_length ; i++)
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				for(unsigned int k = 0 ; k < actual_z_length ; k++)
				{
					data[i][j][k].Imaginary(Randomizer::Uniform(min,max));
				}
			}
		}
	}
	void Signal3D::Increment(const unsigned int& x_index,const unsigned int& y_index,const unsigned int& z_index,const ComplexNumber& number){data[x_index][y_index][z_index] += number;}
	void Signal3D::Decrement(const unsigned int& x_index,const unsigned int& y_index,const unsigned int& z_index,const ComplexNumber& number){data[x_index][y_index][z_index] -= number;}
	void Signal3D::operator*(const double& factor)
	{
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				for(unsigned int k = 0 ; k < actual_z_length ; k++)
				{
					data[i][j][k] *= factor;
				}
			}
		}
	}
	void Signal3D::operator/(const double& factor)
	{
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				for(unsigned int k = 0 ; k < actual_z_length ; k++)
				{
					data[i][j][k] /= factor;
				}
			}
		}
	}
	/*void Signal2D::Transpose()
	{
		ComplexNumber** transposed_data = new ComplexNumber*[y_length];
		for(unsigned int i = 0 ; i < y_length ; i++)
		{
			transposed_data[i] = new ComplexNumber[x_length];
		}
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				transposed_data[j][i] = data[i][j];
			}
			delete [] data[i];
		}
		delete [] data;
		data = transposed_data;
		unsigned int temp = x_length;
		x_length = y_length;
		y_length = temp;
		temp = actual_x_length;
		actual_x_length = actual_y_length;
		actual_y_length = temp;
	}*/
	Signal3D Signal3D::Fourier() const{return FourierAndInverse(false);}
	Signal3D Signal3D::InverseFourier() const{return FourierAndInverse(true);}
	/*unsigned int Signal2D::PowerSpectrum(double*& spectrum) const
	{
		// get the two-dimensional power array
		double** power_array = new double*[x_length];
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			power_array[i] = new double[y_length];
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				power_array[i][j] = data[i][j].SquaredAmplitude();
			}
		}
		// arrange into a one-dimensional power spectrum array
		unsigned int spectrum_size = (unsigned int)ceil(sqrt(x_length*x_length + y_length*y_length));
		spectrum = new double[spectrum_size];
		unsigned int* counts = new unsigned int[spectrum_size];
		for(unsigned int i = 0 ; i < spectrum_size ; i++)
		{
			spectrum[i] = 0.0;
			counts[i] = 0;
		}
		unsigned int index = 0;
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{	
				index = (unsigned int)floor(sqrt(i*i + j*j) + 0.5);
				spectrum[index] += power_array[i][j];
				counts[index]++;
			}
			delete [] power_array[i];
		}
		delete [] power_array;
		// get the average value for each index in the spectrum array
		for(unsigned int i = 0 ; i < spectrum_size ; i++)
		{
			if(counts[i] == 0)			continue;
			spectrum[i] = spectrum[i]/(double)counts[i];
		}
		delete [] counts;
		return spectrum_size;
	}*/
	/*void Signal2D::ShiftCenter()
	{
		ComplexNumber** shifted_data = new ComplexNumber*[x_length];
		unsigned int x_shift = x_length/2;
		unsigned int y_shift = y_length/2;
		unsigned int i_index = 0;
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			i_index = (i + x_shift)%x_length;
			shifted_data[i_index] = new ComplexNumber[y_length];
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				shifted_data[i_index][(j + y_shift)%y_length] = data[i][j];
			}
		}
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			delete [] data[i];
		}
		delete [] data;
		data = shifted_data;
	}
	void Signal2D::AutoCorrelation(Signal2D& autocorrelation,const bool& shift_center) const
	{
		// the Fourier transform of a signal's autocorrelation 
		// is the product of the conjugate of its Fourier transform 
		// by its Fourier transform
		Signal2D transform = Fourier();
		Signal2D product(x_length,y_length);
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				product(i,j,ComplexNumber(transform(i,j).SquaredAmplitude(),0.0));
			}
		}
		// inverse tranform the product to get the autocorrelation
		autocorrelation = product.InverseFourier();
		// shift its center if required
		if(shift_center)		autocorrelation.ShiftCenter();
	}
	double Signal2D::RadialPowerDistribution(const unsigned int& resolution,double*& distribution) const
	{
		// this function assumes that the signal is center shifted
		distribution = new double[resolution];
		unsigned int* counts = new unsigned int[resolution];
		for(unsigned int i = 0 ; i < resolution ; i++)
		{
			distribution[i] = 0.0;
			counts[i] = 0;
		}
		// 0.5 factor because the signal is center shifted and max i/j are at half x/y lengths
		double x = 0.5*(double)x_length;
		double y = 0.5*(double)y_length;
		double delta = sqrt(x*x + y*y)/(double)resolution;
		unsigned int index = 0;
		double x_shift = 0.5*(double)x_length;
		double y_shift = 0.5*(double)y_length;
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{	
				x = double(i) - x_shift;
				y = double(j) - y_shift;
				index = (unsigned int)floor(sqrt(x*x + y*y)/delta);
				distribution[index] += data[i][j].SquaredAmplitude();
				counts[index]++;
			}
		}
		// get the average value for each bin in the distribution array
		for(unsigned int i = 0 ; i < resolution ; i++)
		{
			if(counts[i] == 0)			continue;
			distribution[i] = distribution[i]/(double)counts[i];
		}
		delete [] counts;
		return delta;
	}
	double Signal2D::AngularPowerDistribution(const unsigned int& resolution,double*& distribution) const
	{
		// this function assumes that the signal is center shifted
		distribution = new double[resolution];
		unsigned int* counts = new unsigned int[resolution];
		for(unsigned int i = 0 ; i < resolution ; i++)
		{
			distribution[i] = 0.0;
			counts[i] = 0;
		}
		double delta = 2.0*PI/(double)resolution;
		unsigned int index = 0;
		double x = 0.0;
		double y = 0.0;
		double x_shift = 0.5*(double)x_length;
		double y_shift = 0.5*(double)y_length;
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{	
				// skip the origin
				if((i == x_length/2) && (j == y_length/2))		continue;
				x = double(i) - x_shift;
				y = double(j) - y_shift;
				index = (unsigned int)floor((atan2(y,x) + PI)/delta);
				distribution[index] += data[i][j].SquaredAmplitude();
				counts[index]++;
			}
		}
		// get the average value for each bin in the distribution array
		for(unsigned int i = 0 ; i < resolution ; i++)
		{
			if(counts[i] == 0)			continue;
			distribution[i] = distribution[i]/(double)counts[i];
		}
		delete [] counts;
		return delta;
	}*/
	/*void Complex2DArray::SumAtProcess(const int& process_id)
	{
		unsigned int entry_count = 2*x_length*y_length;
		double* entries = new double[entry_count];
		double* sums = 0;
		if(EZ::Parallel::EZMPI::ProcessID() == process_id)		sums = new double[entry_count];
		unsigned int index = 0;
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			for(unsigned int j = 0 ; j < y_length ; j++)
			{
				entries[index++] = data[i][j].Real();
				entries[index++] = data[i][j].Imaginary();
			}
		}
		EZ::Parallel::EZMPI::Reduce(entries,sums,entry_count,MPI_DOUBLE,MPI_SUM,process_id);
		delete [] entries;
		EZ::Parallel::EZMPI::Barrier();
		if(EZ::Parallel::EZMPI::ProcessID() == process_id)
		{
			index = 0;
			for(unsigned int i = 0 ; i < x_length ; i++)
			{
				for(unsigned int j = 0 ; j < y_length ; j++)
				{
					data[i][j].Real(sums[index++]);
					data[i][j].Imaginary(sums[index++]);
				}
			}
			delete [] sums;
		}
		EZ::Parallel::EZMPI::Barrier();
	}*/
	/*double Signal2D::Test(const unsigned int& x_size,const unsigned int& y_size)
	{
		Signal2D signal(x_size,y_size);
		signal.Randomize(-100.0,100.0);
		//printf("signal:\n");
		//for(unsigned int i = 0 ; i < x_size ; i++)
		//{
			//for(unsigned int j = 0 ; j < y_size ; j++)
			//{
				//printf("\t%f + %f i\n",signal(i,j).Real(),signal(i,j).Imaginary());
			//}
		//}
		//use_fft = false;
		//Signal2D direct_transform = signal.Fourier();
		//use_fft = true;
		Signal2D fft_transform = signal.Fourier();
		//printf("transform:\n");
		//for(unsigned int i = 0 ; i < x_size ; i++)
		//{
			//for(unsigned int j = 0 ; j < y_size ; j++)
			//{
				//printf("\t%f + %f i\n",fft_transform(i).Real(),fft_transform(i).Imaginary());
			//}
		//}
		//use_fft = false;
		//Signal2D direct_direct_recovery = direct_transform.InverseFourier();
		//Signal2D direct_fft_recovery = fft_transform.InverseFourier();
		//use_fft = true;
		//Signal2D fft_direct_recovery = direct_transform.InverseFourier();
		Signal2D fft_fft_recovery = fft_transform.InverseFourier();
		//printf("recovered signal:\n");
		//for(unsigned int i = 0 ; i < x_size ; i++)
		//{
			//for(unsigned int j = 0 ; j < y_size ; j++)
			//{
				//printf("\t%f + %f i : %f + %f i\n",signal(i).Real(),signal(i).Imaginary(),recovery(i).Real(),recovery(i).Imaginary());
			//}
		//}
		ComplexNumber error;
		//double direct_direct_rms_error = 0.0;
		//double fft_direct_rms_error = 0.0;
		//double direct_fft_rms_error = 0.0;
		double fft_fft_rms_error = 0.0;
		for(unsigned int i = 0 ; i < x_size ; i++)
		{
			for(unsigned int j = 0 ; j < y_size ; j++)
			{
				//error = signal(i) - direct_direct_recovery(i);
				//direct_direct_rms_error += error.SquaredAmplitude();
				//error = signal(i) - fft_direct_recovery(i);
				//fft_direct_rms_error += error.SquaredAmplitude();
				//error = signal(i) - direct_fft_recovery(i);
				//direct_fft_rms_error += error.SquaredAmplitude();
				error = signal(i,j) - fft_fft_recovery(i,j);
				fft_fft_rms_error += error.SquaredAmplitude();
			}
		}
		//direct_direct_rms_error = sqrt(direct_direct_rms_error)/(double)size;
		//fft_direct_rms_error = sqrt(fft_direct_rms_error)/(double)size;
		//direct_fft_rms_error = sqrt(direct_fft_rms_error)/(double)size;
		fft_fft_rms_error = sqrt(fft_fft_rms_error/(double)x_size/(double)y_size);
		//printf("direct-direct recovery rms error : %e\n",direct_direct_rms_error);
		//printf("fft-direct recovery rms error : %e\n",fft_direct_rms_error);
		//printf("direct-fft recovery rms error : %e\n",direct_fft_rms_error);
		printf("fft-fft recovery rms error : %e\n",fft_fft_rms_error);
		return fft_fft_rms_error;
	}*/
	/*void Signal2D::WriteVTK(const char* filename,const double& width,const double& height) const
	{
		double dx = width/(double)(actual_x_length - 1);
		double dy = height/(double)(actual_y_length - 1);
		unsigned int point_count = actual_x_length*actual_y_length;
		FILE* file = fopen(filename,"w");
		fprintf(file,"# vtk DataFile Version 3.0\n");
		fprintf(file,"Two-dimensional signal\n");
		fprintf(file,"ASCII\n");
		fprintf(file,"DATASET UNSTRUCTURED_GRID\n");
		fprintf(file,"POINTS %u double\n",point_count);
		for(unsigned int i = 0 ; i < actual_x_length ; i++)
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				fprintf(file,"%15.10f\t%15.10f\t%15.10f\n",i*dx,j*dy,sqrt(data[i][j].SquaredAmplitude()));
			}
		}

		// write the cell connectivity and their types
		unsigned int cell_count = 2*(actual_x_length - 1)*(actual_y_length - 1);
		fprintf(file,"CELLS %u %u\n",cell_count,4*cell_count);
		for(unsigned int i = 0 ; i < (actual_x_length - 1) ; i++)
		{
			for(unsigned int j = 0 ; j < (actual_y_length - 1) ; j++)
			{
				fprintf(file,"3\t%u\t%u\t%u\n",i*actual_y_length + j,(i + 1)*actual_y_length + j,i*actual_y_length + j + 1);
				fprintf(file,"3\t%u\t%u\t%u\n",(i + 1)*actual_y_length + j,(i + 1)*actual_y_length + j + 1,i*actual_y_length + j + 1);
			}
		}
		fprintf(file,"CELL_TYPES %d\n",cell_count);
		for(unsigned int i = 0 ; i < cell_count ; i++)
		{
			fprintf(file,"5\n");
		}
		// write the cell data
		fprintf(file,"POINT_DATA %d\n",point_count);
		fprintf(file,"SCALARS real double\n");
		fprintf(file,"LOOKUP_TABLE default\n");
		for(unsigned int i = 0 ; i < actual_x_length ; i++)
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				fprintf(file,"%e\n",data[i][j].Real());
			}
		}
		fprintf(file,"SCALARS imaginary double\n");
		fprintf(file,"LOOKUP_TABLE default\n");
		for(unsigned int i = 0 ; i < actual_x_length ; i++)
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				fprintf(file,"%e\n",data[i][j].Imaginary());
			}
		}
		fprintf(file,"SCALARS log_amplitude double\n");
		fprintf(file,"LOOKUP_TABLE default\n");
		for(unsigned int i = 0 ; i < actual_x_length ; i++)
		{
			for(unsigned int j = 0 ; j < actual_y_length ; j++)
			{
				fprintf(file,"%e\n",log(sqrt(data[i][j].SquaredAmplitude())));
			}
		}
		fclose(file);
	}*/
	void Signal3D::Initialize()
	{
		actual_x_length = 0;
		actual_y_length = 0;
		actual_z_length = 0;
		x_length = 0;
		y_length = 0;
		z_length = 0;
		data = 0;
	}
	void Signal3D::ClearData()
	{
		if(data == 0)		return;
		for(unsigned int i = 0 ; i < x_length ; i++)
		{
			if(data[i] != 0)
			{
				for(unsigned int j = 0 ; j < y_length ; j++)
				{
					if(data[i][j] != 0)		delete [] data[i][j];
				}
				delete [] data[i];
			}
		}
		delete [] data;
		data = 0;
		x_length = 0;
		y_length = 0;
		z_length = 0;
		actual_x_length = 0;
		actual_y_length = 0;
		actual_z_length = 0;
	}
	/*Signal2D Signal2D::FourierAndInverse(const bool& inverse) const
	{
		Signal2D column_transform(x_length,y_length);
		// get the one dimensional FFT transforms of signal columns
		if(use_fft)
		{
			for(unsigned int i = 0 ; i < x_length ; i++)
			{
				FFT(data[i],column_transform.data[i],y_length,inverse);
			}
		}
		else
		{
			for(unsigned int i = 0 ; i < x_length ; i++)
			{
				DFT(data[i],column_transform.data[i],y_length,inverse);
			}
		}
		// transpose the signal so that we can pass the rows to FFT
		column_transform.Transpose();
		// transform the transposed transform columns (original transform rows)
		Signal2D transform(y_length,x_length);
		if(use_fft)
		{
			for(unsigned int i = 0 ; i < y_length ; i++)
			{
				FFT(column_transform.data[i],transform.data[i],x_length,inverse);
			}
		}
		else
		{
			for(unsigned int i = 0 ; i < x_length ; i++)
			{
				DFT(column_transform.data[i],transform.data[i],x_length,inverse);
			}
		}
		// transpose the transform back
		transform.Transpose();
		transform.actual_x_length = actual_x_length;
		transform.actual_y_length = actual_y_length;
		return transform;
	}*/
}


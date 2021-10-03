// CASDevTest.cpp
// Copyright (c) 2021 Ahmed M. Hussein (amhussein4@gmail.com)
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

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "CASDev.h"
#include "Signal.h"
#include "Image.h"

using namespace CASDev;

void TestFourierTransform1D(const unsigned int& data_size,const unsigned int& test_count)
{
	double error = 0.0;
	double max_error = 0.0;
	for(unsigned int i = 0 ; i < test_count ; i++)
	{
		printf("test %u of %u:\n",i + 1,test_count);
		error = Signal1D::Test(data_size);
		if(error > max_error)		max_error = error;
	}
	printf("max error = %e\n",error);
}

void Test1DSignalAutocorrelation(const unsigned int& data_size)
{
	Signal1D signal(10);
	signal(0,ComplexNumber(1.0,0.0));
	signal(1,ComplexNumber(0.0,0.0));
	signal(2,ComplexNumber(0.0,0.0));
	signal(3,ComplexNumber(1.0,0.0));
	signal(4,ComplexNumber(0.0,0.0));
	signal(5,ComplexNumber(0.0,0.0));
	signal(6,ComplexNumber(0.0,0.0));
	signal(7,ComplexNumber(0.0,0.0));
	signal(8,ComplexNumber(0.0,0.0));
	signal(9,ComplexNumber(0.0,0.0));
	Signal1D autocorelation;
	signal.AutoCorrelation(autocorelation);
	printf("got %u autocorr values\n",autocorelation.Length());
	for(unsigned int i = 0 ; i < autocorelation.Length() ; i++)
	{
		printf("@ %u : %e\t:\t%e\n",i,autocorelation(i).Amplitude(),autocorelation(i).Phase());
	}

	//Signal1D signal(data_size);
	//signal.Randomize(-100.0,100.0);
	//Signal1D autocorelation;
	//signal.AutoCorrelation(autocorelation);
	//printf("got %u autocorr values\n",autocorelation.Length());
	//for(unsigned int i = 0 ; i < autocorelation.Length() ; i++)
	//{
	//	printf("@ %u : %e\t:\t%e\n",i,autocorelation(i).Amplitude(),autocorelation(i).Phase());
	//}
}

void TestFourierTransform2D(const unsigned int& x_data_size,const unsigned int& y_data_size,const unsigned int& test_count)
{
	double error = 0.0;
	double max_error = 0.0;
	for(unsigned int i = 0 ; i < test_count ; i++)
	{
		printf("test %u of %u:\n",i + 1,test_count);
		error = Signal2D::Test(x_data_size,y_data_size);
		if(error > max_error)		max_error = error;
	}
	printf("max error = %e\n",error);
}

void TestPNG(const char* source_name,const char* target_name)
{
	PNG source;
	source.Read(source_name);
	source.Write(target_name);
}

bool ImageFourier(const char* source_name,const char* target_name,const char& channel)
{
	PNG source;
	source.Read(source_name);
	Signal2D source_signal;
	source_signal.ReadFromPNG(source,channel);
	PNG target(source.Width(),source.Height(),1);
	printf("source channel count = %u\n",source.ChannelCount());
	source_signal.WriteToPNG(target,'A',0);
	target.Write(target_name);
	Signal2D transform = source_signal.Fourier();
	PNG transform_png(source.Width(),source.Height(),1);
	transform.WriteToPNG(transform_png,'A',2);
	transform_png.Write("transform.png");
	//transform.WriteVTK("transform.vtk",source.Width(),source.Height());
	/*double* spectrum = 0;
	unsigned int spectrum_size = transform.PowerSpectrum(spectrum);
	for(unsigned int i = 0 ; i < spectrum_size ; i++)
	{
		printf("%u : %e\n",i,spectrum[i]);
	}
	delete [] spectrum;*/
	return true;
}

void ReconstructPNG(const char* source_base_name)
{
	char filename[128];
	sprintf(filename,"%s.png",source_base_name);
	printf("reading %s\n",filename);
	PNG source;
	source.Read(filename);
	Signal2D source_R;
	Signal2D source_G;
	Signal2D source_B;
	Signal2D source_A;
	printf("populating signals\n");
	source_R.ReadFromPNG(source,'R');
	source_G.ReadFromPNG(source,'G');
	source_B.ReadFromPNG(source,'B');
	source_A.ReadFromPNG(source,'A');
	PNG signal_R_png(source.Width(),source.Height(),1);
	PNG signal_G_png(source.Width(),source.Height(),1);
	PNG signal_B_png(source.Width(),source.Height(),1);
	PNG signal_A_png(source.Width(),source.Height(),1);
	printf("writing signals\n");
	source_R.WriteToPNG(signal_R_png,'A',2);
	source_G.WriteToPNG(signal_G_png,'A',2);
	source_B.WriteToPNG(signal_B_png,'A',2);
	source_A.WriteToPNG(signal_A_png,'A',2);
	sprintf(filename,"%s_R_signal.png",source_base_name);
	signal_R_png.Write(filename);
	sprintf(filename,"%s_G_signal.png",source_base_name);
	signal_G_png.Write(filename);
	sprintf(filename,"%s_B_signal.png",source_base_name);
	signal_B_png.Write(filename);
	sprintf(filename,"%s_A_signal.png",source_base_name);
	signal_A_png.Write(filename);
	printf("transforming signals\n");
	Signal2D transform_R = source_R.Fourier();
	Signal2D transform_G = source_G.Fourier();
	Signal2D transform_B = source_B.Fourier();
	Signal2D transform_A = source_A.Fourier();
	printf("writing transforms\n");
	PNG transform_R_png(source.Width(),source.Height(),1);
	PNG transform_G_png(source.Width(),source.Height(),1);
	PNG transform_B_png(source.Width(),source.Height(),1);
	PNG transform_A_png(source.Width(),source.Height(),1);
	transform_R.WriteToPNG(transform_R_png,'A',2);
	transform_G.WriteToPNG(transform_G_png,'A',2);
	transform_B.WriteToPNG(transform_B_png,'A',2);
	transform_A.WriteToPNG(transform_A_png,'A',2);
	sprintf(filename,"%s_R_transform.png",source_base_name);
	transform_R_png.Write(filename);
	sprintf(filename,"%s_G_transform.png",source_base_name);
	transform_G_png.Write(filename);
	sprintf(filename,"%s_B_transform.png",source_base_name);
	transform_B_png.Write(filename);
	sprintf(filename,"%s_A_transform.png",source_base_name);
	transform_A_png.Write(filename);
	printf("inverting transforms\n");
	Signal2D inverse_R = transform_R.InverseFourier();
	Signal2D inverse_G = transform_G.InverseFourier();
	Signal2D inverse_B = transform_B.InverseFourier();
	Signal2D inverse_A = transform_A.InverseFourier();
	printf("writing inverses\n");
	PNG inverse_R_png(source.Width(),source.Height(),1);
	PNG inverse_G_png(source.Width(),source.Height(),1);
	PNG inverse_B_png(source.Width(),source.Height(),1);
	PNG inverse_A_png(source.Width(),source.Height(),1);
	inverse_R.WriteToPNG(inverse_R_png,'A',2);
	inverse_G.WriteToPNG(inverse_G_png,'A',2);
	inverse_B.WriteToPNG(inverse_B_png,'A',2);
	inverse_A.WriteToPNG(inverse_A_png,'A',2);
	sprintf(filename,"%s_R_inverse.png",source_base_name);
	inverse_R_png.Write(filename);
	sprintf(filename,"%s_G_inverse.png",source_base_name);
	inverse_G_png.Write(filename);
	sprintf(filename,"%s_B_inverse.png",source_base_name);
	inverse_B_png.Write(filename);
	sprintf(filename,"%s_A_inverse.png",source_base_name);
	inverse_A_png.Write(filename);
	printf("reconstructing signal\n");
	PNG target(source.Width(),source.Height(),4);
	inverse_R.WriteToPNG(target,'R',2);
	inverse_G.WriteToPNG(target,'G',2);
	inverse_B.WriteToPNG(target,'B',2);
	inverse_A.WriteToPNG(target,'A',2);
	sprintf(filename,"%s_reconstruction.png",source_base_name);
	target.Write(filename);
}

void ImageAutoCorrelation(const char* source_base_name,const unsigned int& resolution)
{
	char filename[128];
	sprintf(filename,"%s.png",source_base_name);
	PNG source;
	source.Read(filename);
	if(source.ChannelCount() == 4)
	{
		Signal2D source_R;
		Signal2D source_G;
		Signal2D source_B;
		Signal2D source_A;
		source_R.ReadFromPNG(source,'R');
		source_G.ReadFromPNG(source,'G');
		source_B.ReadFromPNG(source,'B');
		source_A.ReadFromPNG(source,'A');
		Signal2D autocorrelation_R;
		source_R.AutoCorrelation(autocorrelation_R);
		Signal2D autocorrelation_G;
		source_G.AutoCorrelation(autocorrelation_G);
		Signal2D autocorrelation_B;
		source_B.AutoCorrelation(autocorrelation_B);
		Signal2D autocorrelation_A;
		source_A.AutoCorrelation(autocorrelation_A);
		PNG autocorrelation_R_png(autocorrelation_R.XLength(),autocorrelation_R.YLength(),3);
		PNG autocorrelation_G_png(autocorrelation_G.XLength(),autocorrelation_G.YLength(),3);
		PNG autocorrelation_B_png(autocorrelation_B.XLength(),autocorrelation_B.YLength(),3);
		PNG autocorrelation_A_png(autocorrelation_A.XLength(),autocorrelation_A.YLength(),3);
		autocorrelation_R.WriteLogAmplitudeToPNG(autocorrelation_R_png);
		autocorrelation_G.WriteLogAmplitudeToPNG(autocorrelation_G_png);
		autocorrelation_B.WriteLogAmplitudeToPNG(autocorrelation_B_png);
		autocorrelation_A.WriteLogAmplitudeToPNG(autocorrelation_A_png);
		sprintf(filename,"%s_R_autocorrelation.png",source_base_name);
		autocorrelation_R_png.Write(filename);
		sprintf(filename,"%s_G_autocorrelation.png",source_base_name);
		autocorrelation_G_png.Write(filename);
		sprintf(filename,"%s_B_autocorrelation.png",source_base_name);
		autocorrelation_B_png.Write(filename);
		sprintf(filename,"%s_A_autocorrelation.png",source_base_name);
		autocorrelation_A_png.Write(filename);
	}
	else if(source.ChannelCount() == 1)
	{
		Signal2D source_A;
		source_A.ReadFromPNG(source,'A');
		Signal2D autocorrelation;
		source_A.AutoCorrelation(autocorrelation);
		PNG autocorrelation_png(autocorrelation.XLength(),autocorrelation.YLength(),3);
		autocorrelation.WriteLogAmplitudeToPNG(autocorrelation_png);
		sprintf(filename,"%s_autocorrelation.png",source_base_name);
		autocorrelation_png.Write(filename);
		double* power_r = 0;
		double* power_theta = 0;
		double delta_r = 0.0;
		double delta_theta = 0.0;
		delta_r = autocorrelation.RadialPowerDistribution(resolution,power_r);
		delta_theta = autocorrelation.AngularPowerDistribution(resolution,power_theta);
		for(unsigned int i = 0 ; i < resolution ; i++)
		{
			printf("%u\t%e\t%e\t%e\t%e\n",i,i*delta_r,power_r[i],-PI + i*delta_theta,power_theta[i]);
		}
		delete [] power_r;
		delete [] power_theta;
	}
}

int main(int argc,char** argv)
{
	if(argc < 3)		return 1;
	//TestFourierTransform1D(atoi(argv[1]),atoi(argv[2]));
	//TestFourierTransform2D(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
	//TestPNG(argv[1],argv[2]);
	//ImageFourier(argv[1],argv[2],argv[3][0]);
	//ReconstructPNG(argv[1]);
	//Test1DSignalAutocorrelation(atoi(argv[1]));
	ImageAutoCorrelation(argv[1],atoi(argv[2]));
	return 0;	
}


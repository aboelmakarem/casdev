// Randomizer.h
// Copyright (c) 2009 Ahmed M. Hussein (amhussein4@gmail.com)
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

#ifndef CASDEV_RANDOMIZER_H_
#define CASDEV_RANDOMIZER_H_

namespace CASDev
{
	class Randomizer
	{
	public:
		~Randomizer();
		void Reset();
		static double Uniform();
		static double Uniform(const double& min,const double& max);
		static int UniformInteger(const int& min,const int& max);
		static int Sign();
		static double Normal();
		static double Normal(const double& mean,const double& standard_deviation);
		static double Exponential(const double& mean);
		static void Shuffle(const unsigned int& size,unsigned int* shuffled_array,const unsigned int& passes = 1);
		static void Seed(unsigned int seed);
	
	private:
		Randomizer();
		Randomizer(const Randomizer& randomizer);
		Randomizer& operator=(const Randomizer& randomizer);
		void Initialize();
		static double MTRNG();
		static unsigned int TemperingShiftU(const unsigned int& y);
 		static unsigned int TemperingShiftS(const unsigned int& y);
 		static unsigned int TemperingShiftT(const unsigned int& y);
 		static unsigned int TemperingShiftL(const unsigned int& y);
 		#define MT_N 624
 		#define MT_M 397
 		static unsigned int MT_i;
 		static unsigned int MT_matrix_A;
 		static unsigned int MT_magnitude[2];
 		static unsigned int MT_numbers[MT_N];
 		static unsigned int MT_upper_mask;
 		static unsigned int MT_lower_mask;
 		static unsigned int MT_tempering_mask_b;
 		static unsigned int MT_tempering_mask_c;
	};
}

#endif



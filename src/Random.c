/*	Random.c
	Ahmed M. Hussein (amhussein4@gmail.com)
	04/30/2022

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

#include "Random.h"
#include "math.h"
#include "time.h"

// An implementation for the Mersenne Twister Random Number Generator
#define MT_N 624
#define MT_M 397

unsigned int MT_i = MT_N + 1;
unsigned int MT_magnitude[2] = {0,2567483615};
unsigned int MT_numbers[MT_N];
unsigned int MT_upper_mask = 2147483648;
unsigned int MT_lower_mask = 2147483647;
unsigned int MT_tempering_mask_b = 2636928640;
unsigned int MT_tempering_mask_c = 4022730752;

unsigned int tempering_shift_u(unsigned int y){return (y >> 11);}
unsigned int tempering_shift_s(unsigned int y){return (y << 7);}
unsigned int tempering_shift_t(unsigned int y){return (y << 15);}
unsigned int tempering_shift_l(unsigned int y){return (y >> 18);}

double MT_random_number_generator()
{
	// Use the Mersenne-Twister algorithm to generate random numbers
	unsigned int y = 0;
	if(MT_i >= MT_N)
	{
		unsigned int k = 0;
		if(MT_i == MT_N + 1)
		{
			rand_seed((unsigned int)time(0));
		}
		for(k = 0; k < MT_N - MT_M ; k++)
		{
			y = (MT_numbers[k] & MT_upper_mask)|(MT_numbers[k + 1] & MT_lower_mask);
			MT_numbers[k] = MT_numbers[k + MT_M]^(y >> 1)^MT_magnitude[y & 0x1];
		}
		for(; k < MT_N - 1; k++)
		{
			y = (MT_numbers[k] & MT_upper_mask)|(MT_numbers[k + 1] & MT_lower_mask);
			MT_numbers[k] = MT_numbers[k + (MT_M - MT_N)]^(y >> 1)^MT_magnitude[y & 0x1];
		}
		y = (MT_numbers[MT_N - 1] & MT_upper_mask)|(MT_numbers[0] & MT_lower_mask);
		MT_numbers[MT_N - 1] = MT_numbers[MT_M - 1]^(y >> 1)^MT_magnitude[y & 0x1];
		MT_i = 0;
	}

	y = MT_numbers[MT_i++];
	y = y^tempering_shift_u(y);
	y = y^(tempering_shift_s(y) & MT_tempering_mask_b);
	y = y^(tempering_shift_t(y) & MT_tempering_mask_c);
	y = y^tempering_shift_l(y);
	return ((double)y)/((unsigned int)0xffffffff);
}

// random number generation library functions
double rand_uniform(){return MT_random_number_generator();}
double rand_uniform_interval(double min,double max)
{
	if(max < min)		return 0.0;
	return (min + rand_uniform()*(max - min));
}
int rand_int(int min,int max)
{
	int accepted = 0;
	int result = 0;
	while(!accepted)
	{
		result = (int)floor(rand_uniform_interval(min - 1,max + 1) + 0.5);
		if(result != (min - 1) && result != (max + 1))		accepted = 1;
	}
	return result;
}
int rand_sign()
{
	if(rand_int(0,1) == 0)		return -1;
	return 1;
}
double rand_std_normal()
{
	// Use Box-Muller transform
	double x1 = 0.0;
	double x2 = 0.0;
	double w = 0.0;
	double y1 = 0.0;
	//double y2 = 0.0;
	do
	{
		x1 = 2.0*rand_uniform() - 1.0;
		x2 = 2.0*rand_uniform() - 1.0;
		w = x1*x1 + x2*x2;
	}
	while(w >= 1.0);
	w = sqrt((-2.0*log(w))/w);
	y1 = x1*w;
	//y2 = x2*w;
	return y1;
}
double rand_normal(double mean,double standard_deviation){return (mean + standard_deviation*rand_std_normal());}
double rand_exponential(double mean){return (-mean*log(1.0 - rand_uniform()));}
void rand_shuffle(unsigned int size,unsigned int* shuffled_array,unsigned int passes)
{
	// This function shuffles an array of unsigned integers that runs from 1 to size
	// The function does not allocate the array, the function user is responsible for allocating it
	for(unsigned int i = 0 ; i < size ; i++)
	{
		shuffled_array[i] = i + 1;
	}
	unsigned int index = 0;
	unsigned int temp = 0;
	for(unsigned int i = 0 ; i < passes ; i++)
	{
		for(unsigned int j = 0 ; j < size ; j++)
		{
			index = rand_int(1,size) - 1;
			temp = shuffled_array[j];
			shuffled_array[j] = shuffled_array[index];
			shuffled_array[index] = temp;
		}
	}
}
void rand_seed(unsigned int seed)
{
	MT_numbers[0] = seed;
	for(MT_i = 1 ; MT_i < MT_N ; MT_i++)
	{
		MT_numbers[MT_i] = (69069*MT_numbers[MT_i - 1]);
	}
}
	

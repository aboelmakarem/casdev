// Random.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 04/30/2022

#ifndef RANDOM_H_
#define RANDOM_H_

// random number generation library functions

double rand_uniform();
double rand_uniform_interval(double min,double max);
int rand_int(int min,int max);
int rand_sign();
double rand_std_normal();
double rand_normal(double mean,double standard_deviation);
double rand_exponential(double mean);
void rand_shuffle(unsigned int size,unsigned int* shuffled_array,unsigned int passes);
void rand_seed(unsigned int seed);

#endif


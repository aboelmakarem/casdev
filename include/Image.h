/*	Image.h
	Ahmed M. Hussein (amhussein4@gmail.com)
	05/12/2022

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

#ifndef IMAGE_H_
#define IMAGE_H_

#include "Matrix.h"

typedef struct Image
{
	unsigned int width;
	unsigned int height;
	// color = 1 if the image is a color image (RGB) and 0 if it is grayscale
	// the encoding is the same in both cases but many operations on gray images are faster
	// gray images are encoded such that their pixels R = G = B
	int gray;
	unsigned char* red;
	unsigned char* green;
	unsigned char* blue;
	unsigned char* alpha;
} Image;

// pixel comparators
typedef unsigned char (*PixelComparator)(unsigned char,unsigned char,const unsigned char*);
unsigned char pixcomp_diff(unsigned char a,unsigned char b,const unsigned char* parameters);
unsigned char pixcomp_match_classify(unsigned char a,unsigned char b,const unsigned char* parameters);

Image* create_image();
Image* clone_image(const Image* image);
void reset_image(Image* image);
void destroy_image(Image* image);
int allocate_image(Image* image,unsigned int width,unsigned int height);
int read_png(Image* image,const char* filename);
int write_png(const Image* image,const char* filename);
int write_png_channel(const Image* image,const char* filename,char channel);
unsigned char img_get(const Image* image,unsigned int row,unsigned int column,char channel);
unsigned char img_get_r(const Image* image,unsigned int row,unsigned int column);
unsigned char img_get_g(const Image* image,unsigned int row,unsigned int column);
unsigned char img_get_b(const Image* image,unsigned int row,unsigned int column);
unsigned char img_get_a(const Image* image,unsigned int row,unsigned int column);
void img_set(Image* image,unsigned int row,unsigned int column,char channel,unsigned char value);
void img_set_r(Image* image,unsigned int row,unsigned int column,unsigned char value);
void img_set_g(Image* image,unsigned int row,unsigned int column,unsigned char value);
void img_set_b(Image* image,unsigned int row,unsigned int column,unsigned char value);
void img_set_a(Image* image,unsigned int row,unsigned int column,unsigned char value);
int img_comp(const Image* A,const Image* B,Image* C,unsigned char* parameters,PixelComparator comparator);
int img_histogram(const Image* image,char channel,double* histogram);
int img_threshold(const Image* source,char channel,Image* target);
int img_conv(const Image* source,const Matrix* kernel,Image* convolution);

#endif



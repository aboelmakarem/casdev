/*	Image.cpp
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

#include "Image.h"
#include "png.h"
#include "BLAS.h"
#include "string.h"
#include "stdlib.h"
#include "float.h"
#include "math.h"

// utility functions
int get_channel_index(char channel)
{
	if((channel == 'R') || (channel == 'r'))			return 0;
	if((channel == 'G') || (channel == 'g'))			return 1;
	if((channel == 'B') || (channel == 'b'))			return 2;
	if((channel == 'A') || (channel == 'a'))			return 3;
	return -1;
}

// pixel comparators
unsigned char pixcomp_diff(unsigned char a,unsigned char b,const unsigned char* parameters)
{
	int difference = a - b;
	// the only parameter used here is the tolerance at parameters location 0
	if(abs(difference) < parameters[0])			return 0;
	else if(difference < -parameters[0])		return (difference + 256);
	return difference;
}
unsigned char pixcomp_match_classify(unsigned char a,unsigned char b,const unsigned char* parameters)
{
	int difference = a - b;
	// the parameters used here are:
	// the tolerance at parameters location 0
	// the color of matching pixels at parameters location 1
	// the color of non-matching pixels at parameters location 2
	if(abs(difference) > parameters[0])			return parameters[2];
	return parameters[1];
}

Image* create_image()
{
	Image* image = (Image*)malloc(sizeof(Image));
	image->width = 0;
	image->height = 0;
	// all images are color images by default
	image->gray = 0;
	image->red = 0;
	image->green = 0;
	image->blue = 0;
	image->alpha = 0;
	return image;
}
Image* clone_image(const Image* image)
{
	Image* clone = create_image();
	allocate_image(clone,image->width,image->height);
	unsigned int copy_size = image->height*image->width*sizeof(unsigned char);
	memcpy(clone->red,image->red,copy_size);
	memcpy(clone->green,image->green,copy_size);
	memcpy(clone->blue,image->blue,copy_size);
	memcpy(clone->alpha,image->alpha,copy_size);
	clone->gray = image->gray;
	return clone;
}
void reset_image(Image* image)
{
	if(image->red != 0)				free(image->red);
	if(image->green != 0)			free(image->green);
	if(image->blue != 0)			free(image->blue);
	if(image->alpha != 0)			free(image->alpha);
	image->red = 0;
	image->green = 0;
	image->blue = 0;
	image->alpha = 0;
	image->width = 0;
	image->height = 0;
	image->gray = 0;
}
void destroy_image(Image* image)
{
	if(image != 0)		reset_image(image);
	free(image);
}
int allocate_image(Image* image,unsigned int width,unsigned int height)
{
	if(image == 0)			return 0;
	if(width == 0)			return 0;
	if(height == 0)			return 0;
	unsigned int size = width*height;
	if((image->width*image->height != size))
	{
		reset_image(image);
		// a Image is always 24-bit deep (3 channels, 8 bits per channel)
		// initialize all pixels to zero by default
		image->red = (unsigned char*)calloc(size,sizeof(unsigned char));
		image->green = (unsigned char*)calloc(size,sizeof(unsigned char));
		image->blue = (unsigned char*)calloc(size,sizeof(unsigned char));
		// the alpha channel needs to be set to 0xff (255) by default
		image->alpha = (unsigned char*)malloc(size*sizeof(unsigned char));
		memset(image->alpha,0xff,size*sizeof(unsigned char));
	}
	image->width = width;
	image->height = height;
	return 1;
}
int read_png(Image* image,const char* filename)
{
	FILE* file = fopen(filename,"rb");
	if(file == 0)		return 0;
	// check PNG signature at the beginning of the file
	png_byte signature[8];
	if(fread(signature,sizeof(png_byte),8,file) != 8)		return 0;
	if(png_sig_cmp(signature,0,8) != 0)
	{
		fclose(file);
		return 0;
	}
	// create and initialize the structure to hold PNG data
	png_struct* png_data = png_create_read_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if(png_data == 0)
	{
		fclose(file);
		return 0;
	}
	// create and initialize the structure to hold PNG specifications
	png_info* png_specs = png_create_info_struct(png_data);
	if(png_specs == 0)
	{
		png_destroy_read_struct(&png_data,0,0);
		fclose(file);
		return 0;
	}
	// set the source file
	png_init_io(png_data,file);
	// skip te signature block
	png_set_sig_bytes(png_data,8);
	// read PNG specifications from source file
	png_read_info(png_data,png_specs);
	// extract PNG specifications from spec structure
	png_uint_32 width = png_get_image_width(png_data,png_specs);
	png_uint_32 height = png_get_image_height(png_data,png_specs);
	png_uint_32 color_type = png_get_color_type(png_data,png_specs);
	png_uint_32 bit_depth = png_get_bit_depth(png_data,png_specs);
	// all image pixels are stored as 8-bit deep pixels, if the source PNG is not, strip it 
	// down to 8 bits. This means that upon reading the PNG data, the pixels will be read 
	// as 8-bit pixels
	if(bit_depth == 16)			png_set_strip_16(png_data);
	// all images are stored as 3-channel RGB images. If the source PNG is not RGB PNG, 
	// handle its color type
	// if it is a pallette color, expand it to RGB color
	if(color_type == PNG_COLOR_TYPE_PALETTE)						png_set_palette_to_rgb(png_data);
	// if it is gray but with depth less than 8-bits, expand it to 8-bits
	if((color_type == PNG_COLOR_TYPE_GRAY) && (bit_depth < 8))		png_set_expand_gray_1_2_4_to_8(png_data);
	// if the PNG transparency is not written in its alpha channel, transform the transparency channel to 
	// an alpha channel
	if(png_get_valid(png_data,png_specs,PNG_INFO_tRNS))				png_set_tRNS_to_alpha(png_data);
	// if, after everything, the color type does not have an alpha channel, fill alpha channel with 
	// 100% opacity
	if((color_type == PNG_COLOR_TYPE_RGB) || (color_type == PNG_COLOR_TYPE_GRAY) || 
		(color_type == PNG_COLOR_TYPE_PALETTE))				png_set_filler(png_data,0xff,PNG_FILLER_AFTER);
	// now the PNG will be read as one with an alpha channel, if it is grayscale, set its reading to RGB
	if((color_type == PNG_COLOR_TYPE_GRAY) || (color_type == PNG_COLOR_TYPE_GRAY_ALPHA))	png_set_gray_to_rgb(png_data);
	// all the above transformations were applied to png_data, update png_specs to reflect these transformations
	// when reading
	png_read_update_info(png_data,png_specs);
	// allocate image data
	allocate_image(image,width,height);
	// read image pixel data
	unsigned char** pixels = (unsigned char**)malloc(height*sizeof(unsigned char*));
	for(unsigned int j = 0 ; j < height ; j++)
	{
		pixels[j] = (unsigned char*)malloc(4*width*sizeof(unsigned char));
	}
	png_read_image(png_data,pixels);
	// set pixel data in image arrays
	unsigned int index = 0;
	int is_gray = 1;
	for(unsigned int j = 0 ; j < height ; j++)
	{
		for(unsigned int i = 0 ; i < width ; i++)
		{
			image->red[index] = pixels[j][4*i];
			image->green[index] = pixels[j][4*i + 1];
			image->blue[index] = pixels[j][4*i + 2];
			image->alpha[index] = pixels[j][4*i + 3];
			// check if image is gray
			if(is_gray)
			{
				if(pixels[j][4*i] != pixels[j][4*i + 1])			is_gray = 0;
				else if(pixels[j][4*i] != pixels[j][4*i + 2])		is_gray = 0;
			}
			index++;
		}
		free(pixels[j]);
	}
	free(pixels);
	if(is_gray)				image->gray = 1;
	png_destroy_read_struct(&png_data,&png_specs,0);
	fclose(file);
	return 1;
}
int write_png(const Image* image,const char* filename)
{
	// create and initialize PNG writer
	png_struct* png_writer = png_create_write_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if(png_writer == 0)					return 0;
	png_info* png_specs = png_create_info_struct(png_writer);
	if(png_specs == 0)
	{
		png_destroy_write_struct(&png_writer,0);
		return 0;
	}
	FILE* file = fopen(filename,"wb");
	if(file == 0)
	{
		png_destroy_write_struct(&png_writer,&png_specs);
		return 0;
	}
	png_init_io(png_writer,file);
	//write PNG header, image is always 8-bit deep and has 4 channels
	png_set_IHDR(png_writer,png_specs,image->width,image->height,8,PNG_COLOR_TYPE_RGB_ALPHA,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
	png_write_info(png_writer,png_specs);
	// now write the PNG
	// unpack image pixels from float encoding
	unsigned char** pixels = (unsigned char**)malloc(image->height*sizeof(unsigned char*));
	unsigned int index = 0;
	for(unsigned int j = 0 ; j < image->height ; j++)
	{
		pixels[j] = (unsigned char*)malloc(4*image->width*sizeof(unsigned char));
		for(unsigned int i = 0 ; i < image->width ; i++)
		{
			pixels[j][4*i] = image->red[index];
			pixels[j][4*i + 1] = image->green[index];
			pixels[j][4*i + 2] = image->blue[index];
			pixels[j][4*i + 3] = image->alpha[index];
			index++;
		}
	}
	png_write_image(png_writer,pixels);
	for(unsigned int j = 0 ; j < image->height ; j++)
	{
		free(pixels[j]);
	}
	free(pixels);
	png_write_end(png_writer,png_specs);
	fclose(file);
	png_destroy_write_struct(&png_writer,&png_specs);
	return 1;
}
int write_png_channel(const Image* image,const char* filename,char channel)
{
	// write only 1 channel of the image
	int channel_index = get_channel_index(channel);
	if(channel_index < 0)								return 0;
	
	// create and initialize PNG writer
	png_struct* png_writer = png_create_write_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if(png_writer == 0)					return 0;
	png_info* png_specs = png_create_info_struct(png_writer);
	if(png_specs == 0)
	{
		png_destroy_write_struct(&png_writer,0);
		return 0;
	}
	FILE* file = fopen(filename,"wb");
	if(file == 0)
	{
		png_destroy_write_struct(&png_writer,&png_specs);
		return 0;
	}
	png_init_io(png_writer,file);
	//write PNG header
	png_set_IHDR(png_writer,png_specs,image->width,image->height,8,PNG_COLOR_TYPE_GRAY,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
	png_write_info(png_writer,png_specs);
	// allocate channel pixel data array
	png_byte** pixel_data = (png_byte**)malloc(image->height*sizeof(png_byte*));
	// the bit depth is always 8 and the channel count is always 1
	// populate the pixel data array
	if(image->gray)
	{
		for(unsigned int j = 0 ; j < image->height ; j++)
		{
			pixel_data[j] = (png_byte*)malloc(image->width*sizeof(png_byte));
			for(unsigned int i = 0 ; i < image->width ; i++)
			{
				if(channel_index == 3)					pixel_data[j][i] = img_get_a(image,j,i);
				else   									pixel_data[j][i] = img_get_r(image,j,i);
			}
		}
	}
	else
	{
		for(unsigned int j = 0 ; j < image->height ; j++)
		{
			pixel_data[j] = (png_byte*)malloc(image->width*sizeof(png_byte));
			for(unsigned int i = 0 ; i < image->width ; i++)
			{
				if(channel_index == 0)					pixel_data[j][i] = img_get_r(image,j,i);
				if(channel_index == 1)					pixel_data[j][i] = img_get_g(image,j,i);
				if(channel_index == 2)					pixel_data[j][i] = img_get_b(image,j,i);
				if(channel_index == 3)					pixel_data[j][i] = img_get_a(image,j,i);
			}
		}
	}
	// now write the PNG
	png_write_image(png_writer,pixel_data);
	for(unsigned int j = 0 ; j < image->height ; j++)
	{
		free(pixel_data[j]);
	}
	free(pixel_data);
	png_write_end(png_writer,png_specs);
	fclose(file);
	png_destroy_write_struct(&png_writer,&png_specs);
	return 1;
}
unsigned char img_get(const Image* image,unsigned int row,unsigned int column,char channel)
{
	if(image->gray)				return img_get_r(image,row,column);
	int channel_index = get_channel_index(channel);
	if(channel_index == 0)		return img_get_r(image,row,column);
	if(channel_index == 1)		return img_get_g(image,row,column);
	if(channel_index == 2)		return img_get_b(image,row,column);
	if(channel_index == 3)		return img_get_a(image,row,column);
	return 0;
}
unsigned char img_get_r(const Image* image,unsigned int row,unsigned int column){return image->red[row*image->width + column];}
unsigned char img_get_g(const Image* image,unsigned int row,unsigned int column){return image->green[row*image->width + column];}
unsigned char img_get_b(const Image* image,unsigned int row,unsigned int column){return image->blue[row*image->width + column];}
unsigned char img_get_a(const Image* image,unsigned int row,unsigned int column){return image->alpha[row*image->width + column];}
void img_set(Image* image,unsigned int row,unsigned int column,char channel,unsigned char value)
{
	image->gray = 0;
	int channel_index = get_channel_index(channel);
	if(channel_index == 0)		img_set_r(image,row,column,value);
	if(channel_index == 1)		img_set_g(image,row,column,value);
	if(channel_index == 2)		img_set_b(image,row,column,value);
	if(channel_index == 3)		img_set_a(image,row,column,value);
}
void img_set_r(Image* image,unsigned int row,unsigned int column,unsigned char value){image->red[row*image->width + column] = value;}
void img_set_g(Image* image,unsigned int row,unsigned int column,unsigned char value){image->green[row*image->width + column] = value;}
void img_set_b(Image* image,unsigned int row,unsigned int column,unsigned char value){image->blue[row*image->width + column] = value;}
void img_set_a(Image* image,unsigned int row,unsigned int column,unsigned char value){image->alpha[row*image->width + column] = value;}
int img_comp(const Image* A,const Image* B,Image* C,unsigned char* parameters,PixelComparator comparator)
{
	// compares image A and B and writes comparison result in image C
	// comparison is dictated by the passed comparator function, the comparator decides what is 
	// written in the comparison result
	// check that all images exist and the operands are of the same size
	if(A == 0)			return 0;
	if(B == 0)			return 0;
	if(C == 0)			return 0;
	if(A->width != B->width)		return 0;
	if(A->height != B->height)		return 0;
	allocate_image(C,A->width,A->height);
	unsigned index = 0;
	if(A->gray && B->gray)
	{
		unsigned char temp = 0;
		for(unsigned int j = 0 ; j < A->height ; j++)
		{
			for(unsigned int i = 0 ; i < A->width ; i++)
			{
				temp = comparator(A->red[index],B->red[index],parameters);
				C->red[index] = temp;
				C->green[index] = temp;
				C->blue[index] = temp;
				index++;
			}
		}
	}
	else
	{
		for(unsigned int j = 0 ; j < A->height ; j++)
		{
			for(unsigned int i = 0 ; i < A->width ; i++)
			{				
				C->red[index] = comparator(A->red[index],B->red[index],parameters);
				C->green[index] = comparator(A->green[index],B->green[index],parameters);
				C->blue[index] = comparator(A->blue[index],B->blue[index],parameters);
				index++;
			}
		}
	}
	return 1;
}
int img_histogram(const Image* image,char channel,double* histogram)
{
	// Compute the histogram of the image's channel and store in histogram array
	// histogram array is assumed to have been allocated before calling this 
	// function. It has to hold at least 256 values.
	if(image == 0)				return 0;
	if(histogram == 0)			return 0;
	int channel_index = get_channel_index(channel);
	if(channel_index < 0)								return 0;
	double pixel_count = (double)image->width*(double)image->height;
	memset(histogram,0,256*sizeof(double));
	unsigned int row_start = 0;
	if(channel_index == 0)
	{
		for(unsigned int j = 0 ; j < image->height ; j++)
		{
			row_start = j*image->width;
			for(unsigned int i = 0 ; i < image->width ; i++)
			{
				histogram[image->red[row_start + i]] += 1.0;
			}
		}
	}
	else if(channel_index == 1)
	{
		for(unsigned int j = 0 ; j < image->height ; j++)
		{
			row_start = j*image->width;
			for(unsigned int i = 0 ; i < image->width ; i++)
			{
				histogram[image->green[row_start + i]] += 1.0;
			}
		}
	}
	else if(channel_index == 2)
	{
		for(unsigned int j = 0 ; j < image->height ; j++)
		{
			row_start = j*image->width;
			for(unsigned int i = 0 ; i < image->width ; i++)
			{
				histogram[image->blue[row_start + i]] += 1.0;
			}
		}
	}
	else if(channel_index == 3)
	{
		for(unsigned int j = 0 ; j < image->height ; j++)
		{
			row_start = j*image->width;
			for(unsigned int i = 0 ; i < image->width ; i++)
			{
				histogram[image->alpha[row_start + i]] += 1.0;
			}
		}
	}
	dscal(256,1,histogram,1.0/pixel_count);
	return 1;
}
int img_threshold(const Image* source,char channel,Image* target)
{
	// Compute threshold value and apply it on one channel of the image
	if(source == 0)				return 0;
	if(target == 0)				return 0;
	int channel_index = get_channel_index(channel);
	if(channel_index < 0)								return 0;
	// compute image histogram
	double* histogram = (double*)malloc(256*sizeof(double));
	if(!img_histogram(source,channel,histogram))
	{
		free(histogram);
		return 0;
	}
	// Use Ostu's algorithm and try different tresholding values
	// Otsu's algorithm minimizes the in-class variance of the thresholding classes. The 
	// problem here is that we don't know in advance the threshold value. Hence, go over 
	// all possible threshold values and compute sum of in-class variances and pick the 
	// value that minimizes that number
	double q1 = 0.0;
	double q2 = 0.0;
	double mu1 = 0.0;
	double mu2 = 0.0;
	double sigma1 = 0.0;
	double sigma2 = 0.0;
	double variance_sum = 0.0;
	double min_variance_sum = DBL_MAX;
	unsigned int min_variance_sum_threshold = 0;
	double pixel_count = (double)source->width*(double)source->height;
	double min_q = 1.0e-3/pixel_count;
	for(unsigned int i = 0 ; i < 256 ; i++)
	{
		// compute class probabilities, means and variances given the threshold
		q1 = 0.0;
		mu1 = 0.0;
		sigma1 = 0.0;
		for(unsigned int j = 0 ; j < i ; j++)
		{
			q1 += histogram[j];
			mu1 += ((double)j)*histogram[j];
			sigma1 += ((double)j)*((double)j)*histogram[j];
		}
		if(q1 < min_q)			continue;
		q2 = 0.0;
		mu2 = 0.0;
		sigma2 = 0.0;
		for(unsigned int j = i ; j < 256 ; j++)
		{
			q2 += histogram[j];
			mu2 += j*histogram[j];
			sigma2 += ((double)j)*((double)j)*histogram[j];
		}
		if(q2 < min_q)			continue;
		// scale means by mu1 and mu2 to account for the conditional probabilities (probability of pixel in class 
		// 1 or 2 having a value of i is not historgram[i], it is histogram[i] divided by probability of being in 
		// class 1 or 2)
		mu1 = mu1/q1;
		mu2 = mu2/q2;
		// compute in-class variances and their sum
		sigma1 = sigma1/q1 - mu1*mu1;
		sigma2 = sigma2/q2 - mu2*mu2;
		variance_sum = sigma1 + sigma2;
		if(variance_sum < min_variance_sum)
		{
			min_variance_sum = variance_sum;
			min_variance_sum_threshold = i;
		}
	}
	free(histogram);
	allocate_image(target,source->width,source->height);
	unsigned int index = 0;
	for(unsigned int j = 0 ; j < source->height ; j++)
	{
		index = j*source->width;
		for(unsigned int i = 0 ; i < source->width ; i++)
		{
			if(img_get(source,j,i,channel) > min_variance_sum_threshold)
			{
				target->red[index] = 255;
				target->green[index] = 255;
				target->blue[index] = 255;
			}
			index++;
		}
	}
	target->gray = 1;
	return 1;
}
int img_conv(const Image* source,const Matrix* kernel,Image* convolution)
{
	// Compute the convolution of the given image with the kernel and write the result 
	// into the convolution image
	if(source == 0)										return 0;
	if(kernel == 0)										return 0;
	if(convolution == 0)								return 0;
	// kernel dimensions are assumed to be odd
	int k_shift = kernel->rows/2;
	int l_shift = kernel->columns/2;
	unsigned int min_j = k_shift;
	unsigned int max_j = source->height - k_shift;
	unsigned int min_i = l_shift;
	unsigned int max_i = source->width - l_shift;
	// allocate convolution values matrix and output image
	allocate_image(convolution,source->width,source->height);
	unsigned int row_start = 0;
	double weight = 0.0;
	unsigned int index = 0;
	if(source->gray)
	{
		Matrix* convolution_values = create_matrix(source->height - 2*k_shift,source->width - 2*l_shift);
		double value = 0.0;
		double min_value = DBL_MAX;
		double max_value = -DBL_MAX;
		for(unsigned int j = min_j ; j < max_j ; j++)
		{
			for(unsigned int i = min_i ; i < max_i ; i++)
			{
				value = 0.0;
				for(unsigned int k = 0 ; k < kernel->rows ; k++)
				{
					// set the start of the data row segment that will be weighted by this kernel row 
					row_start = (j + k - k_shift)*source->width + i - l_shift;
					for(unsigned int l = 0 ; l < kernel->columns ; l++)
					{
						weight = mat_get(kernel,k,l);
						value += weight*source->red[row_start + l];
					}
				}
				if(value < min_value)			min_value = value;
				else if(value > max_value)		max_value = value;
				mat_set(convolution_values,j - k_shift,i - l_shift,value);
			}
		}
		// compute bias and scaling factors
		double factor = 255.0/(max_value - min_value);
		// if the convolution values are already within the allowable pixel value range, do not 
		// scale or shift them
		if((min_value >= 0.0) && (max_value <= 255.0))
		{
			factor = 1.0;
			min_value = 0.0;
		}
		// set the convolution image data
		unsigned char conv_value = 0;
		for(unsigned int j = min_j ; j < max_j ; j++)
		{
			row_start = j*source->width;
			for(unsigned int i = min_i ; i < max_i ; i++)
			{
				index = row_start + i;
				conv_value = (unsigned char)floor(factor*(mat_get(convolution_values,j - k_shift,i - l_shift) - min_value) + 0.5);
				convolution->red[index] = conv_value;
				convolution->green[index] = conv_value;
				convolution->blue[index] = conv_value;
			}
		}
		destroy_matrix(convolution_values);
		convolution->gray = 1;
	}
	else
	{
		Matrix* convolution_rvalues = create_matrix(source->height - 2*k_shift,source->width - 2*l_shift);
		Matrix* convolution_gvalues = create_matrix(source->height - 2*k_shift,source->width - 2*l_shift);
		Matrix* convolution_bvalues = create_matrix(source->height - 2*k_shift,source->width - 2*l_shift);
		double rvalue = 0.0;
		double min_rvalue = DBL_MAX;
		double max_rvalue = -DBL_MAX;
		double gvalue = 0.0;
		double min_gvalue = DBL_MAX;
		double max_gvalue = -DBL_MAX;
		double bvalue = 0.0;
		double min_bvalue = DBL_MAX;
		double max_bvalue = -DBL_MAX;
		for(unsigned int j = min_j ; j < max_j ; j++)
		{
			for(unsigned int i = min_i ; i < max_i ; i++)
			{
				rvalue = 0.0;
				gvalue = 0.0;
				bvalue = 0.0;
				for(unsigned int k = 0 ; k < kernel->rows ; k++)
				{
					// set the start of the data row segment that will be weighted by this kernel row 
					row_start = (j + k - k_shift)*source->width + i - l_shift;
					for(unsigned int l = 0 ; l < kernel->columns ; l++)
					{
						weight = mat_get(kernel,k,l);
						index = row_start + l;
						rvalue += weight*source->red[index];
						gvalue += weight*source->green[index];
						bvalue += weight*source->blue[index];
					}
				}
				if(rvalue < min_rvalue)				min_rvalue = rvalue;
				else if(rvalue > max_rvalue)		max_rvalue = rvalue;
				if(gvalue < min_gvalue)				min_gvalue = gvalue;
				else if(gvalue > max_gvalue)		max_gvalue = gvalue;
				if(bvalue < min_bvalue)				min_bvalue = bvalue;
				else if(bvalue > max_bvalue)		max_bvalue = bvalue;
				mat_set(convolution_rvalues,j - k_shift,i - l_shift,rvalue);
				mat_set(convolution_gvalues,j - k_shift,i - l_shift,gvalue);
				mat_set(convolution_bvalues,j - k_shift,i - l_shift,bvalue);
			}
		}
		// compute bias and scaling factors
		double rfactor = 255.0/(max_rvalue - min_rvalue);
		double gfactor = 255.0/(max_gvalue - min_gvalue);
		double bfactor = 255.0/(max_bvalue - min_bvalue);
		// if the convolution values are already within the allowable pixel value range, do not 
		// scale or shift them
		if((min_rvalue >= 0.0) && (max_rvalue <= 255.0))
		{
			rfactor = 1.0;
			min_rvalue = 0.0;
		}
		if((min_gvalue >= 0.0) && (max_gvalue <= 255.0))
		{
			gfactor = 1.0;
			min_gvalue = 0.0;
		}
		if((min_bvalue >= 0.0) && (max_bvalue <= 255.0))
		{
			bfactor = 1.0;
			min_bvalue = 0.0;
		}
		// set the convolution image data
		for(unsigned int j = min_j ; j < max_j ; j++)
		{
			row_start = j*source->width;
			for(unsigned int i = min_i ; i < max_i ; i++)
			{
				index = row_start + i;
				convolution->red[index] = (unsigned char)floor(rfactor*(mat_get(convolution_rvalues,j - k_shift,i - l_shift) - min_rvalue) + 0.5);
				convolution->green[index] = (unsigned char)floor(gfactor*(mat_get(convolution_gvalues,j - k_shift,i - l_shift) - min_gvalue) + 0.5);
				convolution->blue[index] = (unsigned char)floor(bfactor*(mat_get(convolution_bvalues,j - k_shift,i - l_shift) - min_bvalue) + 0.5);
			}
		}
		destroy_matrix(convolution_rvalues);
		destroy_matrix(convolution_gvalues);
		destroy_matrix(convolution_bvalues);
	}
	return 1;
}




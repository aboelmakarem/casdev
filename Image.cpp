// Image.cpp
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

#include "Image.h"
#include "png.h"
#include "string.h"

namespace CASDev
{
	PNG::PNG(){Initialize();}
	PNG::PNG(const PNG& png){*this = png;}
	PNG::PNG(const unsigned int& target_width,const unsigned int& target_height,const unsigned int& target_channel_count)
	{
		Initialize();
		Allocate(target_width,target_height,target_channel_count);
	}
	PNG::~PNG(){Reset();}
	PNG& PNG::operator=(const PNG& png)
	{
		Allocate(png.width,png.height,png.channel_count);
		if(channel_count == 4)
		{
			for(unsigned int i = 0 ; i < height ; i++)
			{
				memcpy(r_channel,png.r_channel,width*sizeof(unsigned char));
				memcpy(g_channel,png.g_channel,width*sizeof(unsigned char));
				memcpy(b_channel,png.b_channel,width*sizeof(unsigned char));
				memcpy(a_channel,png.a_channel,width*sizeof(unsigned char));
			}
		}
		else if(channel_count == 3)
		{
			for(unsigned int i = 0 ; i < height ; i++)
			{
				memcpy(r_channel,png.r_channel,width*sizeof(unsigned char));
				memcpy(g_channel,png.g_channel,width*sizeof(unsigned char));
				memcpy(b_channel,png.b_channel,width*sizeof(unsigned char));
			}
		}
		else if(channel_count == 1)
		{
			for(unsigned int i = 0 ; i < height ; i++)
			{
				memcpy(a_channel,png.a_channel,width*sizeof(unsigned char));
			}
		}
		return *this;
	}
	void PNG::Reset()
	{
		ClearData();
		Initialize();
	}
	void PNG::Allocate(const unsigned int& target_width,const unsigned int& target_height,const unsigned int& target_channel_count)
	{
		Reset();
		channel_count = target_channel_count;
		width = target_width;
		height = target_height;
		if(channel_count == 4)
		{
			r_channel = new unsigned char*[height];
			g_channel = new unsigned char*[height];
			b_channel = new unsigned char*[height];
			a_channel = new unsigned char*[height];
			for(unsigned int i = 0 ; i < height ; i++)
			{
				r_channel[i] = new unsigned char[width];
				g_channel[i] = new unsigned char[width];
				b_channel[i] = new unsigned char[width];
				a_channel[i] = new unsigned char[width];
				memset(r_channel[i],0,width*sizeof(unsigned char));
				memset(g_channel[i],0,width*sizeof(unsigned char));
				memset(b_channel[i],0,width*sizeof(unsigned char));
				memset(a_channel[i],0,width*sizeof(unsigned char));
			}
		}
		else if(channel_count == 3)
		{
			r_channel = new unsigned char*[height];
			g_channel = new unsigned char*[height];
			b_channel = new unsigned char*[height];
			a_channel = 0;
			for(unsigned int i = 0 ; i < height ; i++)
			{
				r_channel[i] = new unsigned char[width];
				g_channel[i] = new unsigned char[width];
				b_channel[i] = new unsigned char[width];
				memset(r_channel[i],0,width*sizeof(unsigned char));
				memset(g_channel[i],0,width*sizeof(unsigned char));
				memset(b_channel[i],0,width*sizeof(unsigned char));
			}
		}
		else if(channel_count == 1)
		{
			r_channel = 0;
			g_channel = 0;
			b_channel = 0;
			a_channel = new unsigned char*[height];
			for(unsigned int i = 0 ; i < height ; i++)
			{
				a_channel[i] = new unsigned char[width];
				memset(a_channel[i],0,width*sizeof(unsigned char));
			}
		}
	}
	bool PNG::Read(const char* filename)
	{
		FILE* source = fopen(filename,"rb");
		if(source == 0)			return false;
		// read the PNG signature
		png_byte signature[8];
		if(fread(signature,sizeof(png_byte),8,source) != 8)
		{
			fclose(source);
			return false;
		}
		if(png_sig_cmp(signature,0,8) != 0)
		{
			fclose(source);
			return false;
		}
		// allocate a struct to read data into
		png_struct* reader = png_create_read_struct(PNG_LIBPNG_VER_STRING,0,0,0);
		if(reader == 0)
		{
			fclose(source);
			return false;
		}
		// allocate information struct
		png_info* information = png_create_info_struct(reader);
		if(information == 0)
		{
			// clear the input struct
			png_destroy_read_struct(&reader,0,0);
			fclose(source);
			return false;
		}
		// initialize the input structure
		png_init_io(reader,source);
		// pass the PNG signature data
		png_set_sig_bytes(reader,8);
		// populate the information struct
		png_read_info(reader,information);
		// read png dimensions, bit depth and channel count
		png_uint_32 target_width = png_get_image_width(reader,information);
		png_uint_32 target_height = png_get_image_height(reader,information);
		png_uint_32 target_bit_depth = png_get_bit_depth(reader,information);
		png_uint_32 target_channel_count = png_get_channels(reader,information);
		// allocate PNG arrays
		Allocate(target_width,target_height,target_channel_count);
		png_byte** rows = new png_byte*[height];
		unsigned int channel_size = target_bit_depth/8;
		unsigned int pixel_size = channel_count*channel_size;
		unsigned int row_size = width*pixel_size;
		for(unsigned int i = 0 ; i < height ; i++)
		{
			rows[i] = new png_byte[row_size];
		}
		// read the image data
		png_read_image(reader,rows);
		unsigned int index = 0;
		for(unsigned int i = 0 ; i < height ; i++)
		{
			if(channel_count == 4)
			{
				for(unsigned int j = 0 ; j < width ; j++)
				{
					index = j*pixel_size;
					memcpy(&r_channel[i][j],&(rows[i][index]),channel_size*sizeof(png_byte));
					memcpy(&g_channel[i][j],&(rows[i][index + channel_size]),channel_size*sizeof(png_byte));
					memcpy(&b_channel[i][j],&(rows[i][index + 2*channel_size]),channel_size*sizeof(png_byte));
					memcpy(&a_channel[i][j],&(rows[i][index + 3*channel_size]),channel_size*sizeof(png_byte));
				}
			}
			else if(channel_count == 3)
			{
				for(unsigned int j = 0 ; j < width ; j++)
				{
					index = j*pixel_size;
					memcpy(&r_channel[i][j],&(rows[i][index]),channel_size*sizeof(png_byte));
					memcpy(&g_channel[i][j],&(rows[i][index + channel_size]),channel_size*sizeof(png_byte));
					memcpy(&b_channel[i][j],&(rows[i][index + 2*channel_size]),channel_size*sizeof(png_byte));
				}
			}
			else if(channel_count == 1)
			{
				for(unsigned int j = 0 ; j < width ; j++)
				{
					index = j*pixel_size;
					memcpy(&a_channel[i][j],&(rows[i][index]),channel_size*sizeof(png_byte));
				}
			}
		}
		// deallocate all allocated arrays and structs
		for(unsigned int i = 0 ; i < height ; i++)
		{
			if(rows[i] != 0)		delete [] rows[i];
		}
		delete [] rows;
		png_destroy_read_struct(&reader,&information,0);
		fclose(source);
		return true;
	}
	bool PNG::Write(const char* filename) const
	{
		// 1 channel data is written to the alpha channel of the output PNG, 3 channel data is written 
		// to the R,G and B channels of the PNG while 4 channel data is written to all 4 PNG channels
		// check the data availability
		if(channel_count == 4)
		{
			if(r_channel == 0)		return false;
			if(g_channel == 0)		return false;
			if(b_channel == 0)		return false;
			if(a_channel == 0)		return false;
		}
		else if(channel_count == 3)
		{
			if(r_channel == 0)		return false;
			if(g_channel == 0)		return false;
			if(b_channel == 0)		return false;
		}
		else if(channel_count == 1)
		{
			if(a_channel == 0)		return false;
		}
		// create a PNG write structure
		png_struct* writer = png_create_write_struct(PNG_LIBPNG_VER_STRING,0,0,0);
		if(writer == 0)			return false;
		png_info* information = png_create_info_struct(writer);
		if(information == 0)
		{
			png_destroy_write_struct(&writer,0);
			return false;
		}
		FILE* target = fopen(filename,"wb");
		if(target == 0)
		{
			png_destroy_write_struct(&writer,0);
			return false;
		}
		// initialize the writer
		png_init_io(writer,target);
		if(channel_count == 4)			png_set_IHDR(writer,information,width,height,8,PNG_COLOR_TYPE_RGB_ALPHA,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
		else if(channel_count == 3)		png_set_IHDR(writer,information,width,height,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
		else if(channel_count == 1)		png_set_IHDR(writer,information,width,height,8,PNG_COLOR_TYPE_GRAY,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
		// write PNG information
		png_write_info(writer,information);
		png_byte** rows = new png_byte*[height];
		for(unsigned int i = 0 ; i < height ; i++)
		{
			rows[i] = new png_byte[width*channel_count];
			if(channel_count == 4)
			{
				for(unsigned int j = 0 ; j < width ; j++)
				{
					rows[i][4*j] = r_channel[i][j];
					rows[i][4*j + 1] = g_channel[i][j];
					rows[i][4*j + 2] = b_channel[i][j];
					rows[i][4*j + 3] = a_channel[i][j];
				}
			}
			else if(channel_count == 3)
			{
				for(unsigned int j = 0 ; j < width ; j++)
				{
					rows[i][3*j] = r_channel[i][j];
					rows[i][3*j + 1] = g_channel[i][j];
					rows[i][3*j + 2] = b_channel[i][j];
				}
			}
			else if(channel_count == 1)
			{
				for(unsigned int j = 0 ; j < width ; j++)
				{
					rows[i][j] = a_channel[i][j];
				}
			}
		}
		png_write_image(writer,rows);
		for(unsigned int i = 0 ; i < height ; i++)
		{
			if(rows[i] != 0)		delete [] rows[i];
		}
		delete [] rows;
		png_write_end(writer,information);
		png_destroy_write_struct(&writer,&information);
		fclose(target);
		return true;
	}
	unsigned int PNG::Width() const{return width;}
	unsigned int PNG::Height() const{return height;}
	unsigned int PNG::ChannelCount() const{return channel_count;}
	unsigned char PNG::R(const unsigned int& row,const unsigned int& column) const{return r_channel[row][column];}
	void PNG::R(const unsigned int& row,const unsigned int& column,const unsigned char& level){r_channel[row][column] = level;}
	unsigned char PNG::G(const unsigned int& row,const unsigned int& column) const{return g_channel[row][column];}
	void PNG::G(const unsigned int& row,const unsigned int& column,const unsigned char& level){g_channel[row][column] = level;}
	unsigned char PNG::B(const unsigned int& row,const unsigned int& column) const{return b_channel[row][column];}
	void PNG::B(const unsigned int& row,const unsigned int& column,const unsigned char& level){b_channel[row][column] = level;}
	unsigned char PNG::A(const unsigned int& row,const unsigned int& column) const{return a_channel[row][column];}
	void PNG::A(const unsigned int& row,const unsigned int& column,const unsigned char& level){a_channel[row][column] = level;}
	void PNG::Initialize()
	{
		width = 0;
		height = 0;
		channel_count = 0;
		r_channel = 0;
		g_channel = 0;
		b_channel = 0;
		a_channel = 0;
	}
	void PNG::ClearData()
	{
		if(channel_count == 4)
		{
			if(r_channel != 0)
			{
				for(unsigned int i = 0 ; i < height ; i++)
				{
					if(r_channel[i] != 0)		delete [] r_channel[i];
				}
				delete [] r_channel;
			}
			if(g_channel != 0)
			{
				for(unsigned int i = 0 ; i < height ; i++)
				{
					if(g_channel[i] != 0)		delete [] g_channel[i];
				}
				delete [] g_channel;
			}
			if(b_channel != 0)
			{
				for(unsigned int i = 0 ; i < height ; i++)
				{
					if(b_channel[i] != 0)		delete [] b_channel[i];
				}
				delete [] b_channel;
			}
			if(a_channel != 0)
			{
				for(unsigned int i = 0 ; i < height ; i++)
				{
					if(a_channel[i] != 0)		delete [] a_channel[i];
				}
				delete [] a_channel;
			}
		}
		else if(channel_count == 3)
		{
			if(r_channel != 0)
			{
				for(unsigned int i = 0 ; i < height ; i++)
				{
					if(r_channel[i] != 0)		delete [] r_channel[i];
				}
				delete [] r_channel;
			}
			if(g_channel != 0)
			{
				for(unsigned int i = 0 ; i < height ; i++)
				{
					if(g_channel[i] != 0)		delete [] g_channel[i];
				}
				delete [] g_channel;
			}
			if(b_channel != 0)
			{
				for(unsigned int i = 0 ; i < height ; i++)
				{
					if(b_channel[i] != 0)		delete [] b_channel[i];
				}
				delete [] b_channel;
			}
		}
		else if(channel_count == 1)
		{
			if(a_channel != 0)
			{
				for(unsigned int i = 0 ; i < height ; i++)
				{
					if(a_channel[i] != 0)		delete [] a_channel[i];
				}
				delete [] a_channel;
			}
		}
		Initialize();
	}
}


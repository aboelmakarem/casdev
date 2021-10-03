// Image.h 
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

#ifndef CASDEV_IMAGE_H_
#define CASDEV_IMAGE_H_

namespace CASDev
{
	class PNG
	{
	public:
		PNG();
		PNG(const PNG& png);
		PNG(const unsigned int& target_width,const unsigned int& target_height,const unsigned int& target_channel_count);
		~PNG();
		PNG& operator=(const PNG& png);
		void Reset();
		void Allocate(const unsigned int& target_width,const unsigned int& target_height,const unsigned int& target_channel_count);
		bool Read(const char* filename);
		bool Write(const char* filename) const;
		unsigned int Width() const;
		unsigned int Height() const;
		unsigned int ChannelCount() const;
		unsigned char R(const unsigned int& row,const unsigned int& column) const;
		void R(const unsigned int& row,const unsigned int& column,const unsigned char& level);
		unsigned char G(const unsigned int& row,const unsigned int& column) const;
		void G(const unsigned int& row,const unsigned int& column,const unsigned char& level);
		unsigned char B(const unsigned int& row,const unsigned int& column) const;
		void B(const unsigned int& row,const unsigned int& column,const unsigned char& level);
		unsigned char A(const unsigned int& row,const unsigned int& column) const;
		void A(const unsigned int& row,const unsigned int& column,const unsigned char& level);

	private:
		void Initialize();
		void ClearData();
		unsigned int channel_count;
		unsigned int width;
		unsigned int height;
		unsigned char** r_channel;
		unsigned char** g_channel;
		unsigned char** b_channel;
		unsigned char** a_channel;
	};
}

#endif


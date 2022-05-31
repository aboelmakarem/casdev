

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "float.h"

extern "C"
{
	#include "Image.h"
}

void TestImage(const char* source_filename,const char* output_basename)
{
	Image* source = create_image();
	read_png(source,source_filename);
	char output_filename[256];
	sprintf(output_filename,"%s_full.png",output_basename);
	write_png(source,output_filename);
	sprintf(output_filename,"%s_red.png",output_basename);
	write_png_channel(source,output_filename,'R');
	sprintf(output_filename,"%s_green.png",output_basename);
	write_png_channel(source,output_filename,'g');
	sprintf(output_filename,"%s_blue.png",output_basename);
	write_png_channel(source,output_filename,'b');
	sprintf(output_filename,"%s_alpha.png",output_basename);
	write_png_channel(source,output_filename,'A');
	destroy_image(source);

	Image* circles = create_image();
	allocate_image(circles,640,480);
	double dx = 0.0;
	double dy = 0.0;
	double r = 0.0;
	for(unsigned int j = 0 ; j < 480 ; j++)
	{
		for(unsigned int i = 0 ; i < 640 ; i++)
		{
			dx = (double)i - 160.0;
			dy = (double)j - 120.0;
			r = sqrt(dx*dx + dy*dy);
			if(r < 50.0)			img_set_r(circles,j,i,255);

			dx = (double)i - 480.0;
			dy = (double)j - 120.0;
			r = sqrt(dx*dx + dy*dy);
			if(r < 50.0)			img_set_g(circles,j,i,255);

			dx = (double)i - 320.0;
			dy = (double)j - 360.0;
			r = sqrt(dx*dx + dy*dy);
			if(r < 50.0)			img_set_b(circles,j,i,255);
		}
	}

	write_png(circles,"circles_full.png");
	write_png_channel(circles,"circles_red.png",'r');
	write_png_channel(circles,"circles_green.png",'G');
	write_png_channel(circles,"circles_blue.png",'B');
	write_png_channel(circles,"circles_alpha.png",'a');
	destroy_image(circles);
}

void TestImageHistogram(const char* source_filename)
{
	Image* source = create_image();
	read_png(source,source_filename);

	double* r_histogram = new double[256];
	double* g_histogram = new double[256];
	double* b_histogram = new double[256];
	img_histogram(source,'r',r_histogram);
	img_histogram(source,'g',g_histogram);
	img_histogram(source,'b',b_histogram);
	destroy_image(source);
	for(unsigned int i = 0 ; i < 256 ; i++)
	{
		printf("%u : %e : %e : %e\n",i,r_histogram[i],g_histogram[i],b_histogram[i]);
	}
	delete [] r_histogram;
	delete [] g_histogram;
	delete [] b_histogram;
}

void TestImageThreshold(const char* source_filename,const char* output_basename,char channel)
{
	Image* source = create_image();
	Image* target = create_image();
	read_png(source,source_filename);
	img_threshold(source,channel,target);
	char output_filename[1024];
	sprintf(output_filename,"%s_threshold.png",output_basename);
	write_png(target,output_filename);
	destroy_image(source);
	destroy_image(target);
}

void TestImageComp(const char* first_filename,const char* second_filename,const char* output_filename)
{
	Image* first = create_image();
	read_png(first,first_filename);
	Image* second = create_image();
	read_png(second,second_filename);
	Image* comparison = create_image();

	unsigned char parameters[3] = {10,0,255};
	img_comp(first,second,comparison,parameters,&pixcomp_match_classify);
	write_png(comparison,output_filename);
	
	destroy_image(first);
	destroy_image(second);
	destroy_image(comparison);
}

void TestConvolution(const char* source_name,const char* target_name)
{
	Image* source = create_image();
	read_png(source,source_name);
	Matrix* kernel = create_matrix(3,3);
	
	mat_set(kernel,0,0,-1.0);
	mat_set(kernel,0,1,-2.0);
	mat_set(kernel,0,2,-1.0);

	mat_set(kernel,1,0,0.0);
	mat_set(kernel,1,1,0.0);
	mat_set(kernel,1,2,0.0);

	mat_set(kernel,2,0,1.0);
	mat_set(kernel,2,1,2.0);
	mat_set(kernel,2,2,1.0);

	Image* convolution = create_image();
	img_conv(source,kernel,convolution);
	write_png(convolution,target_name);

	destroy_image(source);
	destroy_image(convolution);
	destroy_matrix(kernel);
}

int main(int argc,char** argv)
{
	if(argc < 4)
	{
		printf("error: missing run time arguments\n");
		printf("usage:testimage source_image_file_name output_basename channel\n");
		return 1;
	}
	TestImage(argv[1],argv[2]);
	TestImageHistogram(argv[1]);
	TestImageThreshold(argv[1],argv[2],argv[3][0]);
	TestImageComp(argv[1],argv[2],argv[3]);
	TestConvolution(argv[1],argv[2]);
	return 0;
}


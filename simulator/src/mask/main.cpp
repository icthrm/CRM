#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include<unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <mrcimg/mrc2img.h>
#include <mrcimg/img_util.h>

using namespace std;

int main(int argc, char **argv)
{
	struct options opts;
	opts.translate = 0;
	if(GetOpts(argc, argv, &opts) < 0){
		return 0;
	}
	
	util::MrcStack mrcs, mrcnew;
    
	if(!mrcs.Open(opts.input)){
		return 0;
	}
	
	if(mrcs.Y() == 1){
		mrcs.GetHeaderReference().ny = mrcs.Z();
		mrcs.GetHeaderReference().nz = 1;
	}
	
	mrcs.CopyToNewStack(mrcnew);
	
	IplImage* mid = mrcs.GetIplImage(mrcs.Size()/2);
	util::ConvertTo1(mid);
	
	IplImage* mask = cvCreateImage(cvGetSize(mid), mid->depth, mid->nChannels);
	cvFillImage(mask, 0);
	cvCircle(mask, cvPoint(mid->width*.5, mid->height*.5), opts.radius, cvScalar(1), -1);
	
	for(int y = 0; y < mid->height; y++){
		float* ptr = (float*)(mid->imageData+y*mid->widthStep);
		float* mptr = (float*)(mask->imageData+y*mask->widthStep);
		for(int x = 0; x < mid->width; x++){
			*ptr = *ptr* *mptr;
			ptr++;
			mptr++;
		}
	}
	
	if(!opts.translate){
		double avg = 0, dev = 0, min = DBL_MAX, max = DBL_MIN;
		long count = 0;
		
		for(int y = 0; y < mid->height; y++){
			float* ptr = (float*)(mid->imageData+y*mid->widthStep);
			float* mptr = (float*)(mask->imageData+y*mask->widthStep);
			for(int x = 0; x < mid->width; x++){
				if(*mptr){
					if(*ptr < min){
						min = *ptr;
					}
					if(*ptr > max){
						max = *ptr;
					}
					avg += *ptr;
					count++;
					dev += *ptr* *ptr;
				}
				ptr++;
				mptr++;
			}
		}
		
		avg /= count;
		dev = sqrt(dev/count-avg*avg);
		
		
// 		util::GetImgAvgSdvCutoffMinMax(mid, &avg, &dev, &min, &max);
		
		mrcnew.SetHeader(util::MrcStack::MODE_FLOAT, avg-3*dev, avg, avg+3*dev);//min, avg, max);
		mrcnew.SetVolumeSize(mid->width, mid->height, 1);
		mrcnew.WriteHeaderToFile(opts.output);
		mrcnew.AppendIplImageToFile(mid);
	}
	else{
		util::ConvertTo1(mid);
		util::SaveImage(mid, opts.output);
// 		cvSaveImage(opts.output, mid);
	}
	
	cvReleaseImage(&mask);
	cvReleaseImage(&mid);
	
	mrcs.Close();
	mrcnew.Close();
}


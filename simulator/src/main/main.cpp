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
	opts.background = 0;
	if(GetOpts(argc, argv, &opts) < 0){
		return 0;
	}
	
	std::cout<<"Simulation for constraint reconstruction."<<std::endl;
	
	util::MrcStack mrcs, mrcnew;
    
	if(!mrcs.Open(opts.input)){
		return 0;
	}
	
	mrcs.CopyToNewStack(mrcnew);
	
	IplImage* mid = mrcs.GetIplImage(mrcs.Size()/2);
	
	double avg, dev, min, max;    
    util::GetImgAvgSdvCutoffMinMax(mid, &avg, &dev, &min, &max);
	
	IplImage* newimg = cvCreateImage(cvSize(opts.size, opts.size), mid->depth, mid->nChannels);
	
	CvRect rect;
	rect.height = mid->height;
	rect.width = mid->width;
	rect.x = (opts.size-mid->width)*.5;
	rect.y = (opts.size-mid->height)*.5;
	cvSetImageROI(newimg, rect);
	cvCopy(mid, newimg);
	
	cvResetImageROI(newimg);
	cvReleaseImage(&mid);
	
	if(!opts.background){
		cv::Mat noise = cv::Mat(newimg->width, newimg->height, CV_32F);
		cv::RNG rnger(cv::getTickCount());
		rnger.fill(noise, cv::RNG::NORMAL, cv::Scalar::all(avg), cv::Scalar::all(dev));
		IplImage ipltemp = noise;
		IplImage* mask = cvCreateImage(cvSize(opts.size, opts.size), newimg->depth, newimg->nChannels);
		cvFillImage(mask, 1);
		cvCircle(mask, cvPoint(opts.size*.5, opts.size*.5), opts.radius, cvScalar(0), -1);
		for(int y = 0; y < newimg->height; y++){
			float* ptr = (float*)(ipltemp.imageData+y*ipltemp.widthStep);
			float* mptr = (float*)(mask->imageData+y*mask->widthStep);
			for(int x = 0; x < newimg->width; x++){
				*ptr = *ptr* *mptr;
				ptr++;
				mptr++;
			}
		}
		
		cvSmooth(&ipltemp, mask, CV_GAUSSIAN, 0, 0, 0.8);
		// 		cvScaleAdd(&ipltemp, cvScalar(MAX/SDV), cpy, cpy);
		cvAdd(mask, newimg, newimg);
		cvReleaseImage(&mask); 
	}//, cpy); 
	else{
		cv::Mat noise = cv::Mat(newimg->width, newimg->height, CV_32F);
		IplImage ipltemp = noise;
		
		cv::RNG rnger(cv::getTickCount());
		int num = rnger.next()%10+100;
		while(num){
			CvPoint pt;
			pt.x = rnger.next()%opts.size;
			pt.y = rnger.next()%opts.size;
			int rad = rnger.next()%int(opts.size*.5)+5;
			if( sqrt((pt.x-opts.size*.5)*(pt.x-opts.size*.5)+(pt.y-opts.size*.5)*(pt.y-opts.size*.5)) > opts.radius+rad && 
				sqrt((pt.x-opts.size*.5)*(pt.x-opts.size*.5)+(pt.y-opts.size*.5)*(pt.y-opts.size*.5)) < opts.size*.5-rad
			){
				cvCircle(&ipltemp, pt, rad, cvScalar(avg+dev*( (rnger.next()%1000)/1000.0 - .5 ) ), -1);
				num--;
			}
			else{
				continue;
			}
		}
		cvAdd(&ipltemp, newimg, newimg);
	}
	
	if(opts.snr > 0){
		cv::Mat noise = cv::Mat(newimg->width, newimg->height, CV_32F);
		cv::RNG rnger(cv::getTickCount());
		rnger.fill(noise, cv::RNG::NORMAL, cv::Scalar::all(0), cv::Scalar::all(dev/pow10(opts.snr*.1)));
		IplImage ipltemp = noise;
		cvAdd(&ipltemp, newimg, newimg);
	}
	
	util::GetImgAvgSdvCutoffMinMax(newimg, &avg, &dev, &min, &max);
	
	mrcnew.SetHeader(util::MrcStack::MODE_FLOAT, min, avg, max);
	mrcnew.SetVolumeSize(opts.size, 1, opts.size);
	mrcnew.WriteHeaderToFile(opts.output);
	mrcnew.SetVolumeSize(opts.size, opts.size, 1);
	mrcnew.AppendIplImageToFile(newimg);
	
// 	util::ConvertTo1(newimg, true);
// 	util::SaveImage(newimg, "test.pgm");
	
	cvReleaseImage(&newimg);
	mrcs.Close();
	mrcnew.Close();
}


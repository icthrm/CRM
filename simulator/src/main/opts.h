#ifndef OPTS_H__
#define OPTS_H__

#include <iostream>
#include <sstream>
#include <cassert>
extern "C"{
#include <getopt.h>
}
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "util/exception.h"

struct options{
	char input[255];
    int size;
    char output[255];
    int radius;
	int background;
    float snr;
};

void Usage(){
	std::cout<<"[-i INPUT FILENAME]\n"<<std::endl;
	std::cout<<"    MRC file as the interested object in simulation\n"<<std::endl;
	std::cout<<"[-o OUTPUT FILENAME]\n"<<std::endl;
	std::cout<<"    MRC filename for result\n"<<std::endl;
	std::cout<<"[-s SIZE]\n"<<std::endl;
	std::cout<<"    Size of the output mrc (SIZE x SIZE)\n"<<std::endl;
	std::cout<<"[-m MASK RADIUS]\n"<<std::endl;
	std::cout<<"    The radius of the centre mask\n"<<std::endl;
	std::cout<<"[-r SNR VALUE]\n"<<std::endl;
	std::cout<<"    The Signal-to-noise ratio of the additional noise\n"<<std::endl;
	std::cout<<"[-b BACKGROUND KIND]\n"<<std::endl;
	std::cout<<"    0 for random noise stain; 1 for black circles\n"<<std::endl;
}

inline int GetOpts(int argc, char **argv, options* opts_){
	
	static struct option longopts[] ={
		{ "input",    	     required_argument,      NULL,              'i' },
		{ "output",       required_argument,      NULL,              'o' },
        { "size",       required_argument,      NULL,              's' },
		{ "mask_radius",        required_argument,      NULL,              'm' },
        { "snr",        required_argument,      NULL,              'r' },
		{ "background",        required_argument,      NULL,              'b' },
       { NULL,              0,                      NULL,               0  }
    };

	int ch;
	while((ch = getopt_long(argc, argv, "hi:o:s:m:r:b:", longopts, NULL)) != -1){
		switch (ch){

		case '?':
			EX_TRACE("Invalid option '%s'.", argv[optind - 1]);
			return -1;

		case ':':
			EX_TRACE("Missing option argument for '%s'.", argv[optind - 1]);
			return -1;

		case 'h':
			Usage();
			return 0;

		case 'i':{
			std::stringstream iss(optarg);
			iss >> opts_->input;
			if(iss.fail()){
				EX_TRACE("Invalid argument '%s'.", optarg);
				return -1;
			}
		}
		break;
		
		case 'o':{
			std::stringstream iss(optarg);
			iss >> opts_->output;
			if(iss.fail()){
				EX_TRACE("Invalid argument '%s'.", optarg);
			}
		}
		break;

		case 's':{
			std::stringstream iss(optarg);
			iss >> opts_->size;
			if(iss.fail()){
				EX_TRACE("Invalid argument '%s'.", optarg);
			}
		}
		break;

		case 'm':{
			std::stringstream iss(optarg);
			iss >> opts_->radius;
			if(iss.fail()){
				EX_TRACE("Invalid argument '%s'.", optarg);
			}
		}
		break;

		case 'r':{
			std::stringstream iss(optarg);
			iss >> opts_->snr;
			if(iss.fail()){
				EX_TRACE("Invalid argument '%s'.", optarg);
			}
		}
		break;
		
		case 'b':{
			std::stringstream iss(optarg);
			iss >> opts_->background;
			if(iss.fail()){
				EX_TRACE("Invalid argument '%s'.", optarg);
			}
		}
		break;

		case 0:
		break;

		default:
			assert(false);
		} //end switch
	} //end while
	return 1;
}

#endif

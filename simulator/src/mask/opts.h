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
    char output[255];
    int radius;
    bool translate;
};

void Usage(){
	std::cout<<"[-i INPUT FILENAME]\n"<<std::endl;
	std::cout<<"    MRC file\n"<<std::endl;
	std::cout<<"[-o OUTPUT FILENAME]\n"<<std::endl;
	std::cout<<"    MRC file\n"<<std::endl;
	std::cout<<"[-r MASK RADIUS]\n"<<std::endl;
	std::cout<<"    The radius of the centre mask\n"<<std::endl;
	std::cout<<"[-t TRANSLATION]\n"<<std::endl;
	std::cout<<"    translate to tiff file\n"<<std::endl;
}

inline int GetOpts(int argc, char **argv, options* opts_){
	
	static struct option longopts[] ={
		{ "input",    	     required_argument,      NULL,              'i' },
		{ "output",       required_argument,      NULL,              'o' },
        { "translation",       required_argument,      NULL,              't' },
		{ "radius",        required_argument,      NULL,              'r' },
       { NULL,              0,                      NULL,               0  }
    };

	int ch;
	while((ch = getopt_long(argc, argv, "hi:o:r:t", longopts, NULL)) != -1){
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

		case 'r':{
			std::stringstream iss(optarg);
			iss >> opts_->radius;
			if(iss.fail()){
				EX_TRACE("Invalid argument '%s'.", optarg);
			}
		}
		break;
		
		case 't':{
			opts_->translate = true;
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

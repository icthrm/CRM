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
	char dir[1024];
	char output[1024];
	//geometry
	float pitch_angle;
	float zshift;
	int thickness;
	float offset;
	//method
	std::string method;
	int mradius;
	//params for iteration
	int iteration;
	float gamma;
};

inline void Usage(){
	std::cout<<"[-dir DIRECTORY NAME]\n"<<std::endl;
	std::cout<<"    dir that incoude a number of particles\n"<<std::endl;
	std::cout<<"[-o OUTPUT FILENAME]\n"<<std::endl;
	std::cout<<"    MRC filename for result\n"<<std::endl;
	std::cout<<"[-r MASK RADIUS]\n"<<std::endl;
	std::cout<<"    the value of radius for mask\n"<<std::endl;
	std::cout<<"[-g O,P,Z,T]\n"<<std::endl;
	std::cout<<"    Geometry information: offset,pitch_angle,zshift,thickness\n"<<std::endl;
	std::cout<<"[-m METHODS (I,R)]\n"<<std::endl;
	std::cout<<"    SART: SART,iteration_number,relax_parameter\n"<<std::endl;
	std::cout<<"    SIRT: SIRT,iteration_number,relax_parameter\n"<<std::endl;
	std::cout<<"[-h]"<<std::endl;
	std::cout<<"    Help Information\n"<<std::endl;
	std::cout<<"EXAMPLES:\n"<<std::endl;
// 	std::cout<<"-i \n"<<std::endl;
}

inline void PrintOpts(const options& opt){
	std::cout<<"pitch_angle =  "<<opt.pitch_angle<<std::endl;
	std::cout<<"zshift = "<<opt.zshift<<std::endl;
	std::cout<<"thickness = "<<opt.thickness<<std::endl;
	std::cout<<"offset = "<<opt.offset<<std::endl;
	std::cout<<"dir = "<<opt.dir<<std::endl;
	std::cout<<"output = "<<opt.output<<std::endl;
	std::cout<<"method = "<<opt.method<<std::endl;
	std::cout<<"iter = "<<opt.iteration<<std::endl;
	std::cout<<"step = "<<opt.gamma<<std::endl;
}

inline void InitOpts(options* opt){
	opt->pitch_angle = 0;
	opt->zshift = 0;
	opt->thickness = 0;
	opt->offset = 0;
	opt->output[0] = '\0';
}

inline int GetOpts(int argc, char **argv, options* opts_){
	
	static struct option longopts[] ={
       { "help",            no_argument,            NULL,              'h' },
		{ "dir",    	     required_argument,      NULL,              'd' },
		{ "mask_radius",       required_argument,      NULL,            'r' },
		{ "mode",        required_argument,      NULL,              'm' },
		{ "output",       required_argument,      NULL,             'o' },
		{ "geometry",        required_argument,      NULL,              'g' },
       { NULL,              0,                      NULL,               0  }
    };

	int ch;
	while((ch = getopt_long(argc, argv, "hd:r:m:o:g:", longopts, NULL)) != -1){
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

		case 'd':{
			std::stringstream iss(optarg);
			iss >> opts_->dir;
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

		case 'g':{ //offset,xaxistilt,zshift,thickness
			std::stringstream iss(optarg);
			std::string tmp;
			getline(iss, tmp, ',');
			opts_->offset = atof(tmp.c_str());

			getline(iss, tmp, ',');
			opts_->pitch_angle = atof(tmp.c_str());

			getline(iss, tmp, ',');
			opts_->zshift = atof(tmp.c_str());

			getline(iss, tmp);
			opts_->thickness = atoi(tmp.c_str());

			if(iss.fail()){
				EX_TRACE("Invalid argument '%s'.", optarg);
				return -1;
			}
		}
		break;

		case 'm':{
			std::stringstream iss(optarg);
			std::string tmp;
			if ( strcmp(optarg,"BPT") && strcmp(optarg,"RP") ){
				getline(iss, opts_->method, ',');

				if(opts_->method == "SIRT"||opts_->method == "SART"){
					getline(iss, tmp, ',');
					opts_->iteration = atoi(tmp.c_str());
					getline(iss, tmp);
					opts_->gamma = atof(tmp.c_str());
				}
			}
			else{
				getline(iss, opts_->method);
			}
			if(iss.fail()){
				EX_TRACE("Invalid argument '%s'.", optarg);
				return -1;
			}
		}
		break;
		
		case 'r':{
			std::stringstream iss(optarg);
			iss >> opts_->mradius;
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
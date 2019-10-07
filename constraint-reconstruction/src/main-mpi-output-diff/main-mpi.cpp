#include "opts.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "mrcmx/mrcstack.h"
// #include "mrcimg/img_util.h"

#define MILLION 1000000

#define PI_180 0.01745329252f

#ifndef PI
#define     PI  3.14159265358979323846
#endif

#define D2R(__ANGLE__) ((__ANGLE__) * PI_180)


struct Coeff{
	union{
		double p[20];
		struct{
			double a[10];
			double b[10];
		};
	};
	
};

bool ReadAngles(std::vector<float>& angles, const char* name)
{
    std::ifstream in(name);
    if(!in.good()){
        return false;
    }

    while(in.good()){
		float val;
		in>>val;
		if(in.fail()) {
			break;
		}
		
        angles.push_back(val);
    }
    in.close();
    return true;
}

void TranslateAngleToCoefficients(const std::vector<float>& angles, const std::vector<float>& xangles, std::vector<Coeff>& coeffs){
	coeffs.resize(angles.size());
	for(int i = 0; i < angles.size(); i++){
		memset(coeffs[i].p, 0, sizeof(double)*20);
		float beta = D2R(angles[i]);
		float alpha = D2R(xangles[i]);
		
		coeffs[i].a[0] = 0; //
		coeffs[i].a[1] = cos(beta); //x
		coeffs[i].a[2] = sin(alpha)*sin(beta); //y
		coeffs[i].a[3] = -cos(alpha)*sin(beta); //z
		coeffs[i].b[0] = 0; //
		coeffs[i].b[1] = 0; //x
		coeffs[i].b[2] = cos(alpha); //y
		coeffs[i].b[3] = sin(alpha); //z
	}
}

/*solve inverse transfroms defined by Geometry; substitute the inversion into coefficients*/
void DecorateCoefficients(std::vector<Coeff>& coeffs, const Geometry& geo)
{
	double alpha = -D2R(geo.pitch_angle), beta = D2R(geo.offset), t = -geo.zshift;
	double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
	double ca2 = ca*ca, sa2 = sa*sa, cb2 = cb*cb, sb2 = sb*sb;
	
	for(int i = 0; i < coeffs.size(); i++){
		double a[10], b[10];
		memcpy(a, coeffs[i].a, sizeof(double)*10);
		memcpy(b, coeffs[i].b, sizeof(double)*10);
		coeffs[i].a[0] = a[0];
		coeffs[i].a[1] = (a[2]*sa*sb+a[3]*ca*sb+a[1]*cb);//*x
		coeffs[i].a[2] = (a[2]*ca-a[3]*sa);//*y
		coeffs[i].a[3] = (a[3]*ca*cb+a[2]*cb*sa-a[1]*sb);//*z
		coeffs[i].a[4] = (a[4]*ca*cb-a[5]*cb*sa-a[6]*sa2*sb-2*a[9]*ca*sa*sb+a[6]*ca2*sb+2*a[8]*ca*sa*sb);//*x*y
		coeffs[i].a[5] = (a[4]*cb2*sa-a[5]*ca*sb2+a[5]*ca*cb2-2*a[7]*cb*sb-a[4]*sa*sb2+2*a[9]*ca2*cb*sb+2*a[8]*cb*sa2*sb+2*a[6]*ca*cb*sa*sb);//*x*z
		coeffs[i].a[6] = (2*a[8]*ca*cb*sa-a[4]*ca*sb+a[5]*sa*sb+a[6]*ca2*cb-a[6]*cb*sa2-2*a[9]*ca*cb*sa);//*y*z
		coeffs[i].a[7] = (a[7]*cb2+a[5]*ca*cb*sb+a[9]*ca2*sb2+a[8]*sa2*sb2+a[4]*cb*sa*sb+a[6]*ca*sa*sb2);//*x^2
		coeffs[i].a[8] = (a[8]*ca2+a[9]*sa2-a[6]*ca*sa);//*y^2
		coeffs[i].a[9] = (a[7]*sb2+a[9]*ca2*cb2+a[8]*cb2*sa2-a[4]*cb*sa*sb+a[6]*ca*cb2*sa-a[5]*ca*cb*sb);//*z^2
		
		coeffs[i].b[0] = b[0];
		coeffs[i].b[1] = (b[2]*sa*sb+b[3]*ca*sb+b[1]*cb);//*x
		coeffs[i].b[2] = (b[2]*ca-b[3]*sa);//*y
		coeffs[i].b[3] = (b[3]*ca*cb+b[2]*cb*sa-b[1]*sb);//*z
		coeffs[i].b[4] = (b[4]*ca*cb-b[5]*cb*sa-b[6]*sa2*sb-2*b[9]*ca*sa*sb+b[6]*ca2*sb+2*b[8]*ca*sa*sb);//*x*y
		coeffs[i].b[5] = (b[4]*cb2*sa-b[5]*ca*sb2+b[5]*ca*cb2-2*b[7]*cb*sb-b[4]*sa*sb2+2*b[9]*ca2*cb*sb+2*b[8]*cb*sa2*sb+2*b[6]*ca*cb*sa*sb);//*x*z
		coeffs[i].b[6] = (2*b[8]*ca*cb*sa-b[4]*ca*sb+b[5]*sa*sb+b[6]*ca2*cb-b[6]*cb*sa2-2*b[9]*ca*cb*sa);//*y*z
		coeffs[i].b[7] = (b[7]*cb2+b[5]*ca*cb*sb+b[9]*ca2*sb2+b[8]*sa2*sb2+b[4]*cb*sa*sb+b[6]*ca*sa*sb2);//*x^2
		coeffs[i].b[8] = (b[8]*ca2+b[9]*sa2-b[6]*ca*sa);//*y^2
		coeffs[i].b[9] = (b[7]*sb2+b[9]*ca2*cb2+b[8]*cb2*sa2-b[4]*cb*sa*sb+b[6]*ca*cb2*sa-b[5]*ca*cb*sb);//*z^2
	}
	
	//considering z_shift
	for(int i = 0; i < coeffs.size(); i++){
		double a[10], b[10];
		memcpy(a, coeffs[i].a, sizeof(double)*10);
		memcpy(b, coeffs[i].b, sizeof(double)*10);
		
		coeffs[i].a[0] = a[0]+a[3]*t+a[9]*t*t;
		coeffs[i].a[1] = a[1]+a[5]*t;//*x
		coeffs[i].a[2] = a[2]+a[6]*t;//*y
		coeffs[i].a[3] = a[3]+2*a[9]*t;//*z
		
		coeffs[i].b[0] = b[0]+b[3]*t+b[9]*t*t;
		coeffs[i].b[1] = b[1]+b[5]*t;//*x
		coeffs[i].b[2] = b[2]+b[6]*t;//*y
		coeffs[i].b[3] = b[3]+2*b[9]*t;//*z
	}
}

void (*funL)(const Coeff&, double, double, double, double*);

void WarpPosition(const Coeff& coeff, double X, double Y, double Z, double* n){
	n[0] = coeff.a[0]+coeff.a[1]*X+coeff.a[2]*Y+coeff.a[3]*Z+coeff.a[4]*X*Y+coeff.a[5]*X*Z+coeff.a[6]*Y*Z+coeff.a[7]*X*X+coeff.a[8]*Y*Y+coeff.a[9]*Z*Z;
	n[1] = coeff.b[0]+coeff.b[1]*X+coeff.b[2]*Y+coeff.b[3]*Z+coeff.b[4]*X*Y+coeff.b[5]*X*Z+coeff.b[6]*Y*Z+coeff.b[7]*X*X+coeff.b[8]*Y*Y+coeff.b[9]*Z*Z;
}

void LinearPosition(const Coeff& coeff, double X, double Y, double Z, double* n){
	n[0] = coeff.a[0]+coeff.a[1]*X+coeff.a[2]*Y+coeff.a[3]*Z;
	n[1] = coeff.b[0]+coeff.b[1]*X+coeff.b[2]*Y+coeff.b[3]*Z;
}

void ValCoef(const Point3DF& origin, const Point3D& coord, const Coeff& coeff, Weight* wt)
{
	double x, y;

	double X, Y, Z, n[2];
	X = coord.x-origin.x; Y = coord.y-origin.y; Z = coord.z-origin.z;
	
// 	funL(coeff, X, Y, Z, n);
	n[0] = coeff.a[0]+coeff.a[1]*X+coeff.a[2]*Y+coeff.a[3]*Z;
	n[1] = coeff.b[0]+coeff.b[1]*X+coeff.b[2]*Y+coeff.b[3]*Z;
	
	x = n[0]+origin.x; y = n[1]+origin.y;

	wt->x_min = floor(x);
	wt->y_min = floor(y);

	wt->x_min_del = x - wt->x_min;
	wt->y_min_del = y - wt->y_min;
}

static std::ofstream diffput;

void Reproject(const Point3DF& origin, const Volume& vol, const Coeff& coeff, Slice& reproj_val, Slice& reproj_wt){

	Point3D coord; int n;

	for(int z = 0; z < vol.height; z++){
		float* vdrefz = vol.data+z*(size_t)vol.width*vol.length;
		coord.z = z+vol.z;
		
		for(int y = 0; y < vol.length; y++){
			float* vdrefy = vdrefz+y*vol.width;
			coord.y = y+vol.y;
			
			for(int x = 0; x < vol.width; x++){
				float* vdrefx = vdrefy+x;
				coord.x = x+vol.x; Weight wt;
				ValCoef(origin, coord, coeff, &wt);

				if(wt.x_min >= 0 && wt.x_min < vol.width && wt.y_min >= 0 && wt.y_min < vol.length){ //(x_min, y_min)
					n = wt.x_min + wt.y_min * vol.width; //index in reproj
					reproj_val.data[n] += (1-wt.x_min_del) * (1-wt.y_min_del) * (*vdrefx);
					reproj_wt.data[n] += (1-wt.x_min_del) * (1-wt.y_min_del);
				}
				if((wt.x_min+1) >= 0 && (wt.x_min+1) < vol.width && wt.y_min >= 0 && wt.y_min < vol.length){ //(x_min+1, y_min)
					n = wt.x_min+1 + wt.y_min * vol.width; //index in reproj
					reproj_val.data[n] += wt.x_min_del * (1-wt.y_min_del) * (*vdrefx);
					reproj_wt.data[n] += wt.x_min_del * (1-wt.y_min_del);
				}
				if(wt.x_min >= 0 && wt.x_min < vol.width && (wt.y_min+1) >= 0 && (wt.y_min+1) < vol.length){ //(x_min, y_min+1)
					n = wt.x_min + (wt.y_min+1) * vol.width; //index in reproj
					reproj_val.data[n] += (1-wt.x_min_del) * wt.y_min_del * (*vdrefx);
					reproj_wt.data[n] += (1-wt.x_min_del) * wt.y_min_del;
				}
				if((wt.x_min+1) >= 0 && (wt.x_min+1) < vol.width && (wt.y_min+1) >= 0 && (wt.y_min+1) < vol.length){ //(x_min+1, y_min+1)
					n = (wt.x_min+1) + (wt.y_min+1) * vol.width; //index in reproj
					reproj_val.data[n] += wt.x_min_del * wt.y_min_del * (*vdrefx);
					reproj_wt.data[n] += wt.x_min_del * wt.y_min_del;
				}
			}
		}
	}
}

inline void BilinearValue(const Slice& slc, const Weight& wt, float* val, float* vwt)
{
	int n;
	if(wt.x_min >= 0 && wt.x_min < slc.width && wt.y_min >= 0 && wt.y_min < slc.height){ //(x_min, y_min)
		n = wt.x_min + wt.y_min * slc.width;
		*val += (1-wt.x_min_del) * (1-wt.y_min_del) * slc.data[n];
		*vwt += (1-wt.x_min_del) * (1-wt.y_min_del);
	}
	if((wt.x_min+1) >= 0 && (wt.x_min+1) < slc.width && wt.y_min >= 0 && wt.y_min < slc.height){ //(x_min+1, y_min)
		n = wt.x_min+1 + wt.y_min * slc.width;
		*val += wt.x_min_del * (1-wt.y_min_del) * slc.data[n];
		*vwt += wt.x_min_del * (1-wt.y_min_del);
	}
	if(wt.x_min >= 0 && wt.x_min < slc.width && (wt.y_min+1) >= 0 && (wt.y_min+1) < slc.height){ //(x_min, y_min+1)
		n = wt.x_min + (wt.y_min+1) * slc.width;
		*val += (1-wt.x_min_del) * wt.y_min_del * slc.data[n];
		*vwt += (1-wt.x_min_del) * wt.y_min_del;
	}
	if((wt.x_min+1) >= 0 && (wt.x_min+1) < slc.width && (wt.y_min+1) >= 0 && (wt.y_min+1) < slc.height){ //(x_min+1, y_min+1)
		n = wt.x_min+1 + (wt.y_min+1) * slc.width;
		*val += wt.x_min_del * wt.y_min_del * slc.data[n];
		*vwt += wt.x_min_del * wt.y_min_del;
	}
}

void BackProject(const Point3DF& origin, MrcStackM& projs, Volume& vol, Coeff coeffv[])
{
	Slice proj(projs.X(), projs.Y()); 
	Point3D coord;

	memset(vol.data, 0, sizeof(float)*vol.length*vol.width*vol.height);
	
	for(int idx = 0; idx < projs.Z(); idx++){
		printf("BPT begin to read %d projection for %d z-coordinate\n", idx, vol.z);
		projs.ReadSliceZ(idx, proj.data);

		for(int z = 0; z < vol.height; z++){
			float* vdrefz = vol.data+z*(size_t)vol.width*vol.length;
			coord.z = z+vol.z;
			
			for(int y = 0; y < vol.length; y++){
				float* vdrefy = vdrefz+y*vol.width;
				coord.y = y+vol.y;
				
				for(int x = 0; x < vol.width; x++){
					coord.x = x+vol.x; 
					Weight wt; float s = 0, c = 0;
					
					ValCoef(origin, coord, coeffv[idx], &wt);
					BilinearValue(proj, wt, &s, &c);          

					if(c){
						*(vdrefy+x) += (float)(s / c);
					}
				}
			}
		}
	}
} 

void UpdateVolumeByProjDiff(const Point3DF& origin, const Slice& diff, Volume& vol, float gamma, const Coeff& coeff){
	
	Point3D coord;

	for(int z = 0; z < vol.height; z++){
		float* vdrefz = vol.data+z*(size_t)vol.width*vol.length;
		coord.z = z+vol.z;
		
		for(int y = 0; y < vol.length; y++){
			float* vdrefy = vdrefz+y*vol.width;
			coord.y = y+vol.y;
			
			for(int x = 0; x < vol.width; x++){
				coord.x = x+vol.x;
				Weight wt; float s = 0, c = 0;
					
				ValCoef(origin, coord, coeff, &wt);
				BilinearValue(diff, wt, &s, &c);   

				if(c){
					*(vdrefy+x) += (float) (s / c)*gamma;
				}
			}
		}
	}
}

void UpdateWeightsByProjDiff(const Point3DF& origin, const Slice& diff, Volume& values, Volume& weights, const Coeff& coeff){
	
	Point3D coord;

	for(int z = 0; z < values.height; z++){
		float* vdrefz = values.data+z*(size_t)values.width*values.length;
		float* wdrefz = weights.data+z*(size_t)weights.width*weights.length;
		coord.z = z+values.z;
		
		for(int y = 0; y < values.length; y++){
			float* vdrefy = vdrefz+y*values.width;
			float* wdrefy = wdrefz+y*weights.width;
			coord.y = y+values.y;
			
			for(int x = 0; x < values.width; x++){
				coord.x = x+values.x;
				Weight wt; float s = 0, c = 0;
					
				ValCoef(origin, coord, coeff, &wt);
				BilinearValue(diff, wt, &s, &c);   

				*(vdrefy+x) += s;
				*(wdrefy+x) += c;
			}
		}
	}
}

void UpdateVolumeByWeights(Volume& vol, Volume& values, Volume& weights, float gamma){
	size_t pxsize = vol.height*vol.width*vol.length;
	
	for(size_t i = pxsize; i--;){
		if(weights.data[i]){
			vol.data[i] += values.data[i]/weights.data[i]*gamma;
		}
	}
}

void SART(int id, const Point3DF& origin, MrcStackM& projs, std::vector<float*> pslices, Volume& vol, Coeff coeffv[], int iteration, float gamma)
{
	int pxsize = projs.X()*projs.Y();
	Slice reproj_val(projs.X(), projs.Y());	//reprojection value
	Slice reproj_wt(projs.X(), projs.Y());	//reprojection weight
	Slice projection(projs.X(), projs.Y());
	
	for(int i = 0; i < iteration; i++){
		for(int idx = 0; idx < projs.Z(); idx++){
			memset(reproj_val.data, 0, sizeof(float)*pxsize);
			memset(reproj_wt.data, 0, sizeof(float)*pxsize);
			
			Reproject(origin, vol, coeffv[idx], reproj_val, reproj_wt);
		
			MPI_Allreduce(MPI_IN_PLACE, reproj_val.data, pxsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
			MPI_Allreduce(MPI_IN_PLACE, reproj_wt.data, pxsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
			
			printf("SART begin to read %d projection for %d z-coordinate (iteraion %d)\n", idx, vol.z, i);
// 			projs.ReadSliceZ(idx, projection.data);
			memcpy(projection.data, pslices[idx], sizeof(float)*projs.X()*projs.Y());
			
			for(int n = 0; n < pxsize; n++){
				if(reproj_wt.data[n] != 0){
					reproj_val.data[n] /= reproj_wt.data[n];
				}
				reproj_val.data[n] = projection.data[n]-reproj_val.data[n];
			}
			
			if(id == 0){
				double rms2 = 0;
				for(int n = 0; n < pxsize; n++){
					rms2 += reproj_val.data[n]*reproj_val.data[n];
				}
				rms2 /= pxsize;
				diffput<<rms2<<" ";
			}
			
			UpdateVolumeByProjDiff(origin, reproj_val, vol, gamma, coeffv[idx]);
		}
	}
}

void SIRT(const Point3DF& origin, MrcStackM& projs, std::vector<float*> pslices, Volume& vol, Coeff coeffv[], int iteration, float gamma)
{
	int pxsize = projs.X()*projs.Y();
	Slice reproj_val(projs.X(), projs.Y());	//reprojection value
	Slice reproj_wt(projs.X(), projs.Y());	//reprojection weight
	Slice projection(projs.X(), projs.Y());
	
	Volume valvol(vol.x, vol.y, vol.z, vol.length, vol.width, vol.height);
	Volume wtvol(vol.x, vol.y, vol.z, vol.length, vol.width, vol.height);
	
	for(int i = 0; i < iteration; i++){
		memset(valvol.data, 0, sizeof(float)*valvol.length*valvol.width*valvol.height);	
		memset(wtvol.data, 0, sizeof(float)*wtvol.length*wtvol.width*wtvol.height);	
		
		for(int idx = 0; idx < projs.Z(); idx++){
			memset(reproj_val.data, 0, sizeof(float)*pxsize);
			memset(reproj_wt.data, 0, sizeof(float)*pxsize);
			
			Reproject(origin, vol, coeffv[idx], reproj_val, reproj_wt);		//vol is not changed during iteration
		
			MPI_Allreduce(MPI_IN_PLACE, reproj_val.data, pxsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
			MPI_Allreduce(MPI_IN_PLACE, reproj_wt.data, pxsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
			
// 			printf("SIRT begin to read %d projection for %d z-coordinate (iteraion %d)\n", idx, vol.z, i);
// 			projs.ReadSliceZ(idx, projection.data);
			memcpy(projection.data, pslices[idx], sizeof(float)*projs.X()*projs.Y());
			
			for(int n = 0; n < pxsize; n++){
				if(reproj_wt.data[n]){
					reproj_val.data[n] /= reproj_wt.data[n];
				}
				reproj_val.data[n] = projection.data[n]-reproj_val.data[n];
			}
			
			UpdateWeightsByProjDiff(origin, reproj_val, valvol, wtvol, coeffv[idx]);
		}
		
		UpdateVolumeByWeights(vol, valvol, wtvol, gamma);
	}
}

struct SysInfo{
	int id;
	int procs;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int namelen;
};

struct Particle{
	MrcStackM projs, mrcvol;
	std::vector<Coeff> params;
	Geometry geo;
	float* data;
	std::vector<float*> pslices;
	long datasize;
	float gamma;
};


int ATOM(Particle& particle, const char* method, int myid, int procs){
	MrcStackM& projs = particle.projs;
	MrcStackM& mrcvol = particle.mrcvol;
	std::vector<float*>& pslices = particle.pslices;

	std::vector<Coeff>& params = particle.params;
	
	Geometry& geo = particle.geo;
	
	int height;
	int zrem = mrcvol.Z()%procs;
	int volz;   //the start slice of reproject per process
	
	if(myid < zrem){
		height = mrcvol.Z()/procs+1;
		volz = height * myid;
	}
	else{
		height = mrcvol.Z()/procs;
		volz = height * myid+zrem;
	}

	Volume vol(0, 0, volz, mrcvol.Y(), mrcvol.X(), height);//, particle.data+mrcvol.X()*volz);
	memcpy(vol.data, particle.data+vol.z*vol.width, vol.height*vol.width*sizeof(float));
	memset(particle.data, 0, particle.datasize*sizeof(float));

	std::cout<<myid<<": ("<<vol.x<<","<<vol.y<<","<<vol.z<<")"<<"&("<<vol.width<<","<<vol.length<<","<<vol.height<<")"<<std::endl;

	Point3DF origin;
	
	origin.x = mrcvol.X()*.5;
	origin.y = mrcvol.Y()*.5;
	origin.z = mrcvol.Z()*.5;

	if(myid == 0){
		printf("origin.x is %f, origin.y is %f, origin.z is %f\n", origin.x, origin.y, origin.z);
	}

	if(strcmp("SART", method) == 0){
		SART(myid, origin, projs, pslices, vol, &params[0], 1, particle.gamma);
	}
	else{
		SIRT(origin, projs, pslices, vol, &params[0], 1, particle.gamma);
	}
	
	memcpy(particle.data+vol.z*vol.width, vol.data, vol.height*vol.width*sizeof(float));
	MPI_Allreduce(MPI_IN_PLACE, particle.data, particle.datasize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );

	return 0;
}

void GetFilesName(const char* dirname, std::vector<std::string>& filenames)
{
	DIR* dir;
    struct dirent* ptr;
    
    dir = opendir(dirname);
	
	while((ptr = readdir(dir)) != NULL){
		if(strcmp(ptr->d_name, ".") == 0){
			continue;
		}
		if(strcmp(ptr->d_name, "..") == 0){
			continue;
		}
		
        std::string name = ptr->d_name;
		filenames.push_back(name);
	}

    closedir(dir);
}


int main(int argc, char *argv[]){
	SysInfo info;

	MPI_Init(&argc, &argv); //parallel init
	MPI_Comm_rank(MPI_COMM_WORLD, &(info.id));
	MPI_Comm_size(MPI_COMM_WORLD, &(info.procs));
	MPI_Get_processor_name(info.processor_name, &(info.namelen));

	options opts;
    InitOpts(&opts);

    if(GetOpts(argc, argv, &opts) <= 0) {
        EX_TRACE("***WRONG INPUT.\n");
        return -1;
    }
    
    	
	if(info.id == 0){
		diffput.open("diff.txt");
	}
    
    if(info.id == 0){
		PrintOpts(opts);
	}
	
	std::vector<std::string> subdirs;
	GetFilesName(opts.dir, subdirs);
	
	if(!subdirs.size()){
		EX_TRACE("Empty input Dir.\n");
        return -1;
	}
	
	std::vector<Particle> particles(subdirs.size());
	
// 	std::cout<<particles.size()<<std::endl;
	
	for(int i = 0; i < subdirs.size(); i++){
		MrcStackM& projs = particles[i].projs;
		MrcStackM& mrcvol = particles[i].mrcvol;
		std::vector<float*>& pslices = particles[i].pslices;
		if(!projs.ReadFile((std::string(opts.dir)+"/"+subdirs[i]+"/proj.mrc").c_str())){
			printf("File %s cannot access.\n", (std::string(opts.dir)+"/"+subdirs[i]+"/proj.mrc").c_str() );
// 			std::cout<<i<<std::endl;
			return -1;
		}
		
		if(info.id == 0){
			projs.ReadHeader();
		}
		MPI_Bcast(&(projs.header), sizeof(MRCheader), MPI_CHAR, 0, MPI_COMM_WORLD);
		
		pslices.resize(projs.Z());
		for(int idx = 0; idx < projs.Z(); idx++){
			pslices[idx] = new float[projs.X()*projs.Y()];
			projs.ReadSliceZ(idx, pslices[idx]);
		}
		
		projs.Close();
		
		mrcvol.InitializeHeader();
		mrcvol.SetSize(projs.X(), projs.Y(), opts.thickness);

		std::vector<float> angles;
		ReadAngles(angles, (std::string(opts.dir)+"/"+subdirs[i]+"/simu.tlt").c_str() );
		
		std::vector<float> xangles(angles.size(), 0.0);
		
		std::vector<Coeff>& params = particles[i].params;
		TranslateAngleToCoefficients(angles, xangles, params);
		
		Geometry& geo = particles[i].geo;
		geo.offset = opts.offset;
		geo.pitch_angle = opts.pitch_angle;
		geo.zshift = opts.zshift;
		particles[i].gamma = opts.gamma;
		
		DecorateCoefficients(params, geo);
		
		particles[i].datasize = mrcvol.X()*mrcvol.Y()*mrcvol.Z();
		particles[i].data = new float[particles[i].datasize];
	}

	int* mask;
	{int width = particles[0].mrcvol.X(); int height = particles[0].mrcvol.Z();
	
	mask = new int[width*height];
	
	memset(mask, 0, sizeof(int)*width*height);
	
	int cx = width*.5; int cy = height*.5;
	
	for(int x = 0; x < width; x++){
		for(int y = 0; y < height; y++){
			if((x-cx)*(x-cx)+(y-cy)*(y-cy) < opts.mradius*opts.mradius){
				mask[x+y*height] = 1;
			}
		}
	}
	}
	
	{	//sparse Kaczmarz-SART
		
	for(int k = 0; k < opts.iteration+1; k++){
		for(int i = 0; i < subdirs.size(); i++){
			//copy for centre particle
			int preidx = i-1 >= 0 ? i-1:subdirs.size()-1; 
			for(int idx = 0; idx < particles[i].datasize; idx++){
				particles[i].data[idx] = particles[i].data[idx]*(1-mask[idx])+particles[preidx].data[idx]*mask[idx];
			}
			
			ATOM(particles[i], opts.method.c_str(), info.id, info.procs);
		}
		diffput<<std::endl;
	}
	
	}
	diffput.close();
	
// 	{	//sparse Kaczmarz-SART
// 	
// 	cvSmooth();
// 	
// 	IplImage* tmp = cvCreateImage(cvSize(opts.thickness, opts.thickness), IPL_DEPTH_32F, 1);
// 	
// 	for(int k = 0; k < opts.iteration+1; k++){
// 		for(int i = 0; i < subdirs.size(); i++){
// 			//copy for centre particle
// 			int preidx = i-1 >= 0 ? i-1:subdirs.size()-1;
// 			
// 			cvFillImage(tmp, CvScalar(0));
// 			
// 			for(int y = 0; y < tmp->height; y++){
// 				float* ptr = (float*)(tmp->imageData+y*tmp->widthStep);
// 				float* mptr = mask+y*tmp->widthStep;
// 				for(int x = 0; x < tmp->width; x++){
// 					
// 					*ptr = *ptr* *mptr;
// 					ptr++;
// 					mptr++;
// 				}
// 			}
// 			
// 			for(int x = 0; x < tmp->width; x++){
// 				for(int y = 0; y < tmp->height; y++){
// 					if(mask[x+y*opts.thickness]){
// 						
// 					}
// 				}
// 			}
// 			
// 			for(int idx = 0; idx < particles[i].datasize; idx++){
// 				particles[i].data[idx] = particles[i].data[idx]*(1-mask[idx])+particles[preidx].data[idx]*mask[idx];
// 			}
// 			
// 			ATOM(particles[i], opts.method.c_str(), info.id, info.procs);
// 		}
// 	}
// 	
// 	}
	

	
// 	{	//EM-SART
// 	
// 	float* tmp = new float[particles[0].datasize];
// 	
// 	for(int k = 0; k < opts.iteration; k++){
// 		for(int i = 0; i < subdirs.size(); i++){
// 			ATOM(particles[i], info.id, info.procs);
// 		}
// 		
// 		memset(tmp, 0, sizeof(float)*particles[0].datasize);
// 		
// 		for(int i = 0; i < subdirs.size(); i++){
// 			//copy for centre particle
// 			for(int idx = 0; idx < particles[i].datasize; idx++){
// 				tmp[idx] += particles[i].data[idx]*mask[idx];
// 			}
// 		}
// 		
// 		for(int idx = 0; idx < particles[0].datasize; idx++){
// 			tmp[idx] /= subdirs.size();
// 		}
// 		
// 		for(int i = 0; i < subdirs.size(); i++){
// 			for(int idx = 0; idx < particles[i].datasize; idx++){
// 				particles[i].data[idx] = particles[i].data[idx]*(1-mask[idx])+tmp[idx]*mask[idx];
// 			}
// 		}
// 	}
// 	delete [] tmp;
// 	
// 	}
	
	MrcStackM& mrcvol = particles.back().mrcvol;//particles[9].mrcvol;//particles.back().mrcvol;
// 	float* data = particles[0].data; //particles.back().data;//particles[9].data; //
	float* rslt = new float[particles[0].datasize];
	memset(rslt, 0, sizeof(int)*particles[0].datasize);
	
	for(int i = 0; i < subdirs.size(); i++){
		for(int idx = 0; idx < particles[0].datasize; idx++){
			rslt[idx] += particles[i].data[idx];
		}
	}
	
	for(int idx = 0; idx < particles[0].datasize; idx++){
		rslt[idx] /= subdirs.size();
	}
	
	float* data = rslt;
	
	mrcvol.WriteToFile(opts.output);
	
	if(info.id == 0){
		mrcvol.WriteHeader();
		Volume vol(0, 0, 0, mrcvol.Y(), mrcvol.X(), mrcvol.Z(), data);
		mrcvol.WriteBlock(vol.z, vol.z+vol.height, 'z', vol.data);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	if(info.id == 0){
		mrcvol.UpdateHeader();
	}

// 	projs.Close();
// 	mrcvol.Close();
	delete [] rslt;
	delete [] mask;
	for(int i = 0; i < particles.size(); i++){
		delete [] particles[i].data;
		for(int j = 0; j < particles[i].projs.Z(); j++){
			delete particles[i].pslices[j];
		}
	}
	
	//delete pslices
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();                 //parallel finish

    return 0;
}


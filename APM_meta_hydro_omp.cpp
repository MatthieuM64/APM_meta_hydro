/*C++ CODE - MANGEAT MATTHIEU - 2025*/
/*ACTIVE POTTS MODEL - HYDRODYNAMICS*/

//////////////////////
///// LIBRAIRIES /////
//////////////////////

//Public librairies.
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>

using namespace std;

//Personal libraries.
#include "lib/special_functions.cpp"

//////////////////////////////
///// MESHGRID FUNCTIONS /////
//////////////////////////////

//Initial condition.
void meshInit(vector<vector<double> > &RHO1, vector<vector<double> > &RHO2, vector<vector<double> > &RHO3, vector<vector<double> > &RHO4, const int &init, const double &rho0, const double &beta, const int &LX, const int &LY, const double &Rd, const double &rhod)
{
	const double mag=0.8; //mag of the liquid (will update quickly to the equilibrium value).
	
	//droplet inside a liquid.
	#pragma omp parallel for default(shared)
	for (int x=0; x<LX; x++)
	{
		for (int y=0; y<LY; y++)
		{
			if (square(x-0.5*LX)+square(y-0.5*LY)<square(Rd)) //inside the droplet.
			{
				RHO1[x][y]=0;
				RHO2[x][y]=rhod*init;
				RHO3[x][y]=rhod*(1-init);
				RHO4[x][y]=0;
			}
			else //outside the droplet.
			{
				RHO1[x][y]=0.25*rho0*(1+3*mag);
				RHO2[x][y]=0.25*rho0*(1-mag);
				RHO3[x][y]=0.25*rho0*(1-mag);
				RHO4[x][y]=0.25*rho0*(1-mag);
			}
		}
	}
}

//Flipping terms.
double Iflip(const double &rj, const double ri, const double &rho, const double &beta)
{
	double mag=4*beta*(rj-ri)/rho;
	return exp(1./(2*rho))*((rj+ri-1./(4*beta))*sinh(mag)-(rj-ri)*cosh(mag));
}

//Finite difference algorithm.
void finiteDiff(vector<vector<double> > &RHO1, vector<vector<double> > &RHO2, vector<vector<double> > &RHO3, vector<vector<double> > &RHO4, const double &Dpara, const double &Dperp, const double &v0, const double &gamma0, const double &beta, const int &LX, const int &LY)
{
	static vector<double> DRHO1(LX*LY,0), DRHO2(LX*LY,0), DRHO3(LX*LY,0), DRHO4(LX*LY,0);
	
	#pragma omp parallel default(shared) 
	{
	
	#pragma omp for schedule(static)
	for (int x=0; x<LX; x++)
	{
		const int xm=x-1+(x==0)*LX, xp=x+1-(x==LX-1)*LX; //x periodic neighbours
		#pragma omp simd
		for (int y=0; y<LY; y++)
		{
			const int ym=y-1+(y==0)*LY, yp=y+1-(y==LY-1)*LY; //y periodic neighbours
			const double rho=RHO1[x][y]+RHO2[x][y]+RHO3[x][y]+RHO4[x][y];
			
			//Flipping terms.
			const double I12=Iflip(RHO1[x][y],RHO2[x][y],rho,beta);
			const double I13=Iflip(RHO1[x][y],RHO3[x][y],rho,beta);
			const double I14=Iflip(RHO1[x][y],RHO4[x][y],rho,beta);
			const double I23=Iflip(RHO2[x][y],RHO3[x][y],rho,beta);
			const double I24=Iflip(RHO2[x][y],RHO4[x][y],rho,beta);
			const double I34=Iflip(RHO3[x][y],RHO4[x][y],rho,beta);
			
			//Hydrodynamic equation (density variations).
			DRHO1[x*LY+y]=Dpara*(RHO1[xm][y]+RHO1[xp][y]-2*RHO1[x][y])+Dperp*(RHO1[x][ym]+RHO1[x][yp]-2*RHO1[x][y]) - v0*(RHO1[xp][y]-RHO1[xm][y]) + gamma0*(I12+I13+I14);
			DRHO2[x*LY+y]=Dperp*(RHO2[xm][y]+RHO2[xp][y]-2*RHO2[x][y])+Dpara*(RHO2[x][ym]+RHO2[x][yp]-2*RHO2[x][y]) - v0*(RHO2[x][yp]-RHO2[x][ym]) + gamma0*(-I12+I23+I24);
			DRHO3[x*LY+y]=Dpara*(RHO3[xm][y]+RHO3[xp][y]-2*RHO3[x][y])+Dperp*(RHO3[x][ym]+RHO3[x][yp]-2*RHO3[x][y]) + v0*(RHO3[xp][y]-RHO3[xm][y]) + gamma0*(-I13-I23+I34);
			DRHO4[x*LY+y]=Dperp*(RHO4[xm][y]+RHO4[xp][y]-2*RHO4[x][y])+Dpara*(RHO4[x][ym]+RHO4[x][yp]-2*RHO4[x][y]) + v0*(RHO4[x][yp]-RHO4[x][ym]) + gamma0*(-I14-I24-I34);
		}
	}
	
	#pragma omp for schedule(static)
	for (int x=0; x<LX; x++)
	{
		#pragma omp simd
		for (int y=0; y<LY; y++)
		{
			//Implement the density variations.
			RHO1[x][y]+=DRHO1[x*LY+y];
			RHO2[x][y]+=DRHO2[x*LY+y];
			RHO3[x][y]+=DRHO3[x*LY+y];
			RHO4[x][y]+=DRHO4[x*LY+y];
		}
	}
	
	}
}

///////////////////////////
///// BASIC FUNCTIONS /////
///////////////////////////

//Total average on the all space.
double average(const vector< vector<double> > &RHO, const int &LX, const int &LY)
{
	double rhoAv=0;
	for (int x=0; x<LX; x++)
	{
		for (int y=0; y<LY; y++)
		{
			rhoAv+=RHO[x][y];
		}
	}
	return double(rhoAv)/(LX*LY);
}

//Export magnetization.
void exportMag(const vector< vector<double> > &RHO1, const vector< vector<double> > &RHO2, const vector< vector<double> > &RHO3, const vector< vector<double> > &RHO4, const double &beta, const double &epsilon, const double &rho0, const int &LX, const int &LY, const double &Rd, const double &rhod, const int &NX, const int &NY, const int &init, const double &t)
{
	const int a=NX/LX;
	
	//Selected magnetization to export.
	static vector<float> mag(LX*LY,0.);
	for (int y=0;y<LY;y++)
	{
		for (int x=0; x<LX; x++)
		{
			double rho=RHO1[a*x][a*y]+RHO2[a*x][a*y]+RHO3[a*x][a*y]+RHO4[a*x][a*y];
			if (init==0)
			{
				mag[y*LX+x]=(4*RHO3[a*x][a*y]-rho)/(3*rho);
			}
			else
			{
				mag[y*LX+x]=(4*RHO2[a*x][a*y]-rho)/(3*rho);
			}
		}
	}
	
	//Write in a binary file.
	int returnSystem=system("mkdir -p data_APM4_dynamics2d/");
	
	stringstream ssMAG;
	ssMAG << "./data_APM4_dynamics2d/APM4_mag_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_Rd=" << Rd << "_rhod=" << rhod << "_init=" << init << "_t=" << t << ".bin";
	string nameMAG = ssMAG.str();
	ofstream fileMAG(nameMAG.c_str(),ios::binary);
	fileMAG.write(reinterpret_cast<const char*>(mag.data()), mag.size()*sizeof(float));
	fileMAG.close();		
}

//Read parameters from command line.
void ReadCommandLine(int argc, char** argv, double &beta, double &rho0, double &epsilon, int &LX, int &LY, double &Rd, double &rhod, int &init, double &dx, double &dt, double &tmax, int &THREAD_NUM)
{
 	for( int i = 1; i<argc; i++ )
	{
		if (strstr(argv[i], "-beta=" ))
		{
			beta=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-rho0=" ))
		{
			rho0=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-epsilon=" ))
		{
			epsilon=atof(argv[i]+9);
		}
		else if (strstr(argv[i], "-LX=" ))
		{
			LX=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-LY=" ))
		{
			LY=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-Rd=" ))
		{
			Rd=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-rhod=" ))
		{
			rhod=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-init=" ))
		{
			init=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-dt=" ))
		{
			dt=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-dx=" ))
		{
			dx=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-tmax=" ))
		{
			tmax=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-threads=" ))
		{
			THREAD_NUM=atoi(argv[i]+9);
		}
		else
		{
			cerr << "BAD ARGUMENT : " << argv[i] << endl;
			abort();
		}
	}
}

/////////////////////
///// MAIN CODE /////
/////////////////////

int main(int argc, char *argv[])
{
	//Physical parameters: beta=inverse temperature, epsilon=self-propulsion parameter, rho0=average density, Rd=droplet radius, rhod=droplet density, LX*LY=size of the box.
	double beta=0.75, epsilon=2.4, rho0=3., Rd=10, rhod=15;
	int LX=200, LY=200;
	
	//Numerical parameters: dt=discrete time, tmax=maximal time, dx=meshgrid size, init=initial condition, THREAD_NUM=number of threads.
	double dt=0.00025, tmax=10000., dx=0.25;
	int init=0, THREAD_NUM=4;

	//Read imported parameters in command line.
	ReadCommandLine(argc,argv,beta,rho0,epsilon,LX,LY,Rd,rhod,init,dx,dt,tmax,THREAD_NUM);
	
	//OpenMP.
	omp_set_dynamic(0);
	omp_set_num_threads(THREAD_NUM);
	const int OMP_MAX_THREADS=omp_get_max_threads();
	cout << OMP_MAX_THREADS << " maximum threads on this node. " << THREAD_NUM << " threads will be used." << endl;
	
	//Grid size.
	const int NX=int(LX/dx), NY=int(LY/dx);
	
	//Global variables.
	const double Dpara=(1+epsilon/3.)*dt/square(dx), Dperp=(1-epsilon/3.)*dt/square(dx), v0=2*epsilon*dt/(3*dx), gamma0=dt;

	//Density and magnetization at each vertex.
	vector<vector<double> > RHO1(NX,vector<double>(NY,0)), RHO2(NX,vector<double>(NY,0)), RHO3(NX,vector<double>(NY,0)), RHO4(NX,vector<double>(NY,0));
	
	//Creation of the file to export global averages.
	const int dossier=system("mkdir -p ./data_APM4_averages/");
	stringstream ssAverages;
	ssAverages << "./data_APM4_averages/APM4_averages_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << ".txt";
	string nameAverages = ssAverages.str();
	ofstream fileAverages(nameAverages.c_str(),ios::trunc);
	fileAverages.precision(6);
	
	//Increments.
	const int Nsteps=int(tmax/dt), DeltaT=int(1./dt);
	
	//Initial condition.
	meshInit(RHO1,RHO2,RHO3,RHO4,init,rho0,beta,NX,NY,Rd/dx,rhod);
	
	//Time evolution.
	for(int t=0;t<=Nsteps;t++)
	{
		//Export data.
		if (t%DeltaT==0 or t==Nsteps)
		{
			exportMag(RHO1,RHO2,RHO3,RHO4,beta,epsilon,rho0,LX,LY,Rd,rhod,NX,NY,init,t*dt);
			
			const double n1=average(RHO1,NX,NY), n2=average(RHO2,NX,NY), n3=average(RHO3,NX,NY), n4=average(RHO4,NX,NY);
			cout << "time=" << t*dt << " -rho=" << n1+n2+n3+n4 << " -n1=" << n1 << " -n2=" << n2 << " -n3=" << n3 << " -n4=" << n4 << running_time.TimeRun(" ") << endl;
			fileAverages <<  t*dt << " " << n1+n2+n3+n4 << " " << n1 << " " << n2 << " " << n3 << " " << n4 << endl;
		}
		
		//At each time-step update RHO and MAG arrays.
		finiteDiff(RHO1,RHO2,RHO3,RHO4,Dpara,Dperp,v0,gamma0,beta,NX,NY);
	}
	return 0;
}

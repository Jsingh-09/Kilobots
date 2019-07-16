/******************************************************************************************************\

This prgram written and maintained by : Katie Copenhagen at UC Merced, Physics.

	This simulates a group of agents that are connected to neighboring 
	   agents by springs which can be fairly long range as long as no
	   other agents are interrupting the space inbetween the agents 
	   connected by the spring.
	
	The agents self propel and can align with eachother, they also have
	   polydispersity and angular noise with a gaussian distribution.

\******************************************************************************************************/

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <time.h>
#include <vector>

using namespace std;

#define N 300		//Number of agents. [0-10000]
#define AA 0		//Alignment interaction magnitude. [0-10]
#define V 0.05		//Agent velocity for self propulsion. [0-0.1]
#define K 0.0		//Spring interaction strength. [0-0.1]
#define NOISE 0.0	//Noise (variance of gaussian distribution for random angular noise). [0-2*pi]
#define RS 0.0		//Variance in agent sizes, for gaussian distribution. [0-RA (defined below)]
#define GAP 0		//The amount of overlap to break springs. [0-1]
#define epsi 0.5      // interaction stength 
#define sigm 1      // diameter of the sphere
#define NU 4       // replusive term power -> for NU <= 3, no thermodynamically stable phases are found. NU = 12 DOES NOT HAVE ANY OUTPUT

#define L 1		//Size of the circle which the agents are initiallized within. [sqrt(N) ish]
#define RA 1		//Average size of the agents, sets lengthscale. [1]
#define T 10000		//Number of timesteps the simulation is run for. [100-100000]
#define SAMPLES 1	//Number of samples run. Must be 1 to output files for making videos. [1-100+]

#define _USE_MATH_DEFINES

#include "particle.cpp"		//Particle class.
#include "clustering.cpp"	//Groups agents into clusters.
#include "opramcl.cpp"		//Calculate the properties of the agnet clusters.

int main(){
	//Seed the random number generator.
	long int letseed=time(NULL);
	srand (letseed);

	float t1=clock();	//Save starting time to calculate the program runtime with.

	FILE *f00;	//For making videos.
	FILE *f01;	//Store parameters.
	f01 = fopen("params.txt", "w");
	FILE *f02;	//Cluster data at simulation finish.
	f02=fopen("data/clusters.dat","w");
	FILE *f03;	//Average order vs. t.
	f03=fopen("data/ordervst.dat","w");
	FILE *f04;	//Agent data at simulation finish.
	f04=fopen("data/agents.dat","w");

	vector<particle> agents(N, particle());		//Vector for storing all the agents.
	for (int s=0; s<SAMPLES; s++){			//Loop through the samples.
		//Initialize agents, random position and velocity direction.
		for (int i=0;i<N;i++){
			agents[i].init();
		}

		for (int t=0; t<T; t++){		//Loop through time.

			//Once enough time has passed to stablize the cluster, spread the sizes, introducing the polydispersity gradually over 100 timesteps.
			if (t>(int)3*(sqrtf((float)N)/((float)(2.0*K))*RA)&&t<(int)3*(sqrtf(N)/((float)(2.0*K))*RA)+100){
				for (int i=0;i<N;i++)
					agents[i].spreadSizes();
					printf("I am here\n");
			}

			//Create neighbor list for all agents.
			for (int i=0; i<N; i++){
				vector<int> temp;	//Temporarily stores the neighbor list.
				temp.reserve(N);	//Pre reserve space to speed up program.

				for (int j=0; j<N; j++){	//Loop through all other particles to check if they should be added to neighborlist.
					float dx=agents[j].getX()-agents[i].getX();
					float dy=agents[j].getY()-agents[i].getY();
					float d=sqrtf(dx*dx+dy*dy);
					if (d<(agents[i].getSize()+agents[j].getSize())+RA*4){
						temp.push_back(j);  //If i and j are close enough together add j to the neighborlist, (within 4 average agent diameters of eachother).
					}
				}
				agents[i].setNeigh(temp); //Set the neighborlist to temp.
				temp.clear();
				agents[i].setOldAng();  //Also save their angles for use in interactions.
			}


			//Loop throgh every particle to calculate interactions.
			for (int i=0; i<N; i++){
				float Ax=0, Ay=0, A, Sx=0, Sy=0, newx, newy, dx, dy, d, newang;	//For storing all the interactions, to calculate new directions for agent i.
				for (int j=0; j<agents[i].getNeigh().size(); j++){		//Loop through j: all of agent i's neighbors.
					//Calculate separation vector between agent i and i's j'th neighbor to use for spring.
					float dx=agents[agents[i].getNeigh()[j]].getX()-agents[i].getX();
					float dy=agents[agents[i].getNeigh()[j]].getY()-agents[i].getY();
					float d=sqrtf(dx*dx+dy*dy);
					//Make sure it doesnt blow up when they are too close.
					if (d<0.001){
						d=0.001;
					}
					if (d>20.0){
						d = 20.0;
					}
					//Spring interaction.
					int kspringtemp=1;	//kspringtemp is for storing whether the spring is interrupted by other cells or not.
					float dx2, dy2, d2;	//For storning separation vector to other neighbors than j.
					if (i!=j){		
						float avgr=(agents[i].getSize()+agents[agents[i].getNeigh()[j]].getSize());
						for (int k=0; k<agents[i].getNeigh().size(); k++){	//Loop through all other neighbors of i that aren't j, they shallt be called k. 
							if (j!=k && i!=k){				//Only do the thing if they are 3 different agents.
								//Separation vector between agent i and it's k'th neighbor.
								dx2=agents[agents[i].getNeigh()[k]].getX()-agents[i].getX();	
								dy2=agents[agents[i].getNeigh()[k]].getY()-agents[i].getY();
								d2=sqrtf(dx2*dx2+dy2*dy2);
								//b basically records the size of the space between agent i and the j'th neighbor, at the relative distance to the k'th neighbor.
								float b=agents[i].getSize()+(agents[agents[i].getNeigh()[j]].getSize()-agents[i].getSize())/(d*d)*(dx2*dx+dy2*dy);
								//Then check if the k'th neighbor is within the GAP parameter fraction of the space between agent i and the j'th neighbor. If it is then interrupt the spring by setting kspringtemp to 0.
								// if (sqrtf(d2*d2-powf(dx*dx2/d+dy*dy2/d,2))<agents[agents[i].getNeigh()[k]].getSize()+GAP*b&&(dx*dx2+dy*dy2)>0){
								// 	kspringtemp=0;
								// }
							}
						}
						//Alignment interaction.
						Ax+=cosf(agents[agents[i].getNeigh()[j]].getOldAng());
						Ay+=sinf(agents[agents[i].getNeigh()[j]].getOldAng());
						// Ax += 0;
						// Ay += 0;
						if (d<avgr){	//Always have the repulsive part of the spring interaction.
							kspringtemp=1;
						}
						//Spring interaction.
						// Sx+=-kspringtemp*(avgr-d)*(dx/d);
						// Sy+=-kspringtemp*(avgr-d)*(dy/d);
						Sx+=kspringtemp*epsi*pow(sigm,NU)*pow((avgr-d),-(NU+1))*(dx/d); // -> cos(angle) = dx/d
						Sy+=kspringtemp*epsi*pow(sigm,NU)*pow((avgr-d),-(NU+1))*(dy/d); // -> sin(angle) = dy/d 
					}
				}

				//Parameters u and v just for calculating the gaussion random variable for noise.
				float u=rand()/(float)RAND_MAX;
				float v=rand()/(float)RAND_MAX;
				//New direction of propulsion for agent i is a gaussian random variable with variance NOISE, plus the angle that comes from the alignment interaction, and it's own travel direction.
				agents[i].setAng(fmodf(NOISE*sqrtf(-2*log(u))*cosf(2*M_PI*v)+atan2(AA*Ay+sinf(agents[i].getAng()),AA*Ax+cosf(agents[i].getAng())),2*M_PI));
				//Set the velocity of the agent to the propulsion with magnitude V in the propulsion direction stored in Ang, and the spring interaction with magnitude K.
				agents[i].setVx(V*cosf(agents[i].getAng())+Sx);
				agents[i].setVy(V*sinf(agents[i].getAng())+Sy);

			}

			//Move agents according to their new directions.
			for (int i=0; i<N; i++)
				agents[i].updatePos();

			int vidtime=10;		//How many timesteps to skip between each frame for the video.
			if (SAMPLES==1 && t%vidtime==0){
				char buff[32];		//buff saves the name of the file to store this frame data in.
				snprintf(buff, sizeof(char)*32, "video/frames/fr%06d", t);
				f00=fopen(buff, "w");	//Files for saving video frame data to make videos (vid.py).
				//Loop through each agent to put data into video frame file.
				for (int i=0;i<N;i++){
					//Find the relative velocity of each agent with all of it's neighbors (will set the color of agents in the video).
					float relv=0;
					for (int j=0;j<agents[i].getNeigh().size();j++){
						//Find how far agent i has moved since the last video frame.
						float mxi=agents[i].getX()-agents[i].getOldX();
						float myi=agents[i].getY()-agents[i].getOldY();
						//Find how far the j'th neighbor has moved since last video frame.
						float mxj=agents[agents[i].getNeigh()[j]].getX()-agents[agents[i].getNeigh()[j]].getOldX();
						float myj=agents[agents[i].getNeigh()[j]].getY()-agents[agents[i].getNeigh()[j]].getOldY();
						relv+=sqrtf(powf((mxj-mxi),2)+powf((myj-myi),2));	//Relative movement of agent i and all neighbors.
					}
					relv=relv/(agents[i].getNeigh().size()*V*vidtime);	//Normalize relative velocity..
					//Write agent data into file for video frame structure 4 columns: relative velocity	X position	Y position	agent size.
					fprintf(f00, "%f	%f	%f	%f\n",relv,agents[i].getX(),agents[i].getY(),agents[i].getSize());
				}
				//Save the current agent positions to use for calculating relative velocities of the next frame.
				for (int i=0; i<N; i++)
					agents[i].storePos();
				fclose(f00);

			}
		}	//End of time loop.

		//Calculate properties of each cluster.
		int labels[N], clusternum;
		clusternum = clustering(agents, labels);		//Function to assign a cluster label to each agent, and returns the total number of clusters.
		for (int i=0; i<N; i++)
			agents[i].setLabel(labels[i]);

		float Order[N],MedianR[N],AvgR[N],StdR[N];		//Arrays for storing the properties of each cluster, only clusternum elements will be saved into these arrays.
		int Size[N];
		opramcl(agents, Size, Order, MedianR, AvgR, StdR, clusternum);	//Calculate and store the size, order, median, mean, and standard deviation of the size of agents within each cluster.

		for (int i=0; i<clusternum; i++){	//Loop through the clusters and store the data about them.
			fprintf(f02,"%d	%f	%f	%f	%f\n", Size[i], Order[i], MedianR[i], AvgR[i], StdR[i]);
		}
		for (int i=0; i<N; i++){	//Save the individual position, direction, size, cluster, and neighbor lists for each agent. 
			fprintf(f04,"%f	%f	%f	%f	%d	%d	",agents[i].getX(),agents[i].getY(),agents[i].getAng(),agents[i].getSize(),agents[i].getLabel(),(int)agents[i].getNeigh().size());
			for (int j=0; j<agents[i].getNeigh().size(); j++)
				fprintf(f04,"%d	",agents[i].getNeigh()[j]);
			fprintf(f04,"\n");
		}
	}

	float t2=clock();
	//Write all of the system parameters and run time to params.txt file.
	fprintf(f01, "N = 	%d \nv = 	%f \nra = 	%f\nrs =	%f\na =	%f\nk =	%f\nnoise = 	%f\nsamples = 	%d\nmaxtime =	%d\nseed = 	%ld \n runtime = %f", N, V, (float)RA, (float)RS, (float)AA, (float)K, (float)NOISE, SAMPLES, T, letseed, (t2-t1)/(float)CLOCKS_PER_SEC);
	fclose(f01);
	return 0;
}

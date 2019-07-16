#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>

void opramcl(vector<particle> agents, int Size[], float Order[], float MedianR[], float AvgR[], float StdR[], int clusternum){
	float txvel[N], tyvel[N], avgr[N], stdr[N];		//For storing the total velocity, average and standard deviation of size to calculate the order parameters.
	int i, size[N];						//Store the size (number of particles in) cluster.

	vector< vector<float> > sizesformed(clusternum); //Sizes of the agents in each cluster for calculating median of the agent size for each cluster.

	//Initialize things to zero for however many clusters there are.
	for (i=0; i<clusternum; i++){
		txvel[i]=0;
		tyvel[i]=0;
		size[i]=0;
		avgr[i]=0;
		stdr[i]=0;
	}

	for (i=0; i<N; i++){	//agents[i].getLable() returns the index for the cluster that agent i is in.
		//Calculate total velocities for finding order parameter.
		txvel[agents[i].getLabel()]+=cosf(agents[i].getAng());
		tyvel[agents[i].getLabel()]+=sinf(agents[i].getAng());
		//Total number of particles in the cluster.
		size[agents[i].getLabel()]++;
		//Average size of the particles in the cluster, and a list of all the sizes for finding the median.
		avgr[agents[i].getLabel()]+=agents[i].getSize();
		sizesformed[agents[i].getLabel()].push_back(agents[i].getSize());
	}

	for (i=0; i<N; i++){
		//Standard deviation of particle size.
		stdr[agents[i].getLabel()]+=powf(agents[i].getSize()-avgr[agents[i].getLabel()]/(float)size[agents[i].getLabel()],2.0);
	}

	//Loop through all the clusters to finalize calculations of parameters and store them in the arrays to be accessed in main.cpp.
	for (i=0; i<clusternum; i++){
		if (txvel[i]*txvel[i]+tyvel[i]*tyvel[i]>0)
			Order[i]=(sqrtf(txvel[i]*txvel[i]+tyvel[i]*tyvel[i])/(float)size[i]);		//Order paramter is magnitude of the sum of the velocities divided by the size of the cluster.
		Size[i]=(size[i]);		//The size is the number of particles in the cluster.

		size_t n = sizesformed[i].size();		//Last index of the cluster, so that when the vector of particle sizes for this cluster is sorted, the n/2 element will be the median.
		sort(sizesformed[i].begin(), sizesformed[i].end());
		//If there are an even number of particles in the cluster then median is the after of the middle 2.
		if (n%2==0)
			MedianR[i]=float(sizesformed[i][n/2 -1] + sizesformed[i][n/2])/2.0;
		else		//If there are an odd number then it is the middle particle.
			MedianR[i]=sizesformed[i][n/2];

		AvgR[i]=(avgr[i]/(float)size[i]);		//Average size of particles in the cluster.
		StdR[i]=(sqrtf(stdr[i]/(float)size[i]));	//Standard deviation of particle sizes in the cluster.
	}
	return;
}

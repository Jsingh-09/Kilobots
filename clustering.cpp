#include <iostream>

void addNewNeighbor(int i, int currentLabel, vector<particle> agents, int labels[]){
	int j;
	if (labels[i]==-1){		//If particle i isnt in a cluster yet add it to the current cluster.
		labels[i]=(currentLabel);	
		for (j=0; j<agents[i].getNeigh().size(); j++){		//Go through all of particle i's neighbors and add them to the same cluster as particle i.
			if (agents[i].getNeigh()[j]!=i)			//Also add all of particle j's neighbors to the cluster, do this recursively until all of the particles in this cluster have a label that is not -1 anymore.
				addNewNeighbor(agents[i].getNeigh()[j], currentLabel, agents, labels);
		}
	}	
	return;		//Once all of the neighbors and neighbors neighbors etc... are added to the current cluster return to start the next cluster.
}

int clustering(vector<particle> agents, int labels[]){		//Assign a label for each particle by cluster.
	int j, currentLabel=0;					//Start with label 0.
	for (j=0; j<N; j++)	
		labels[j]=-1;					//Set all the labels for all particles to -1 to begin with.

	for (j=0; j<N; j++){					//Loops through all the particles to assign them to a cluster
		if (labels[j]==-1){				//If particle j is not yet in a cluster, add it to the current cluster.
			addNewNeighbor(j, currentLabel, agents, labels);		//Add particle j and it's neighbors to the current cluster.
			currentLabel++;				//Increase the current cluster label by 1 to start labeling the particles in the next cluster in the next loop through.
		}
	}
	return currentLabel;		//currentLabel at the end is the total number of cluster, so return that.
}

#include <iostream>

using namespace std;
class particle {
	private:
		//Position of particle.
		float X;
		float OldX;
		float Y;
		float OldY;
		//Velocity of particle.
		float Vx;
		float Vy;
		//Direction of self propulsion of particle.
		float Ang;
		float OldAng;
		//Current size of particle, and equilibrium size after spread has happened.
		float Size;
		float Sizeeq;
		//Lable of the cluster that particle is in.
		int Label;
		//Neighborlist.
		vector<int> Neigh;
		
	public:
		//Set the variables.
		void setAng(float ang){Ang=ang;}
		void setOldAng(void){OldAng=Ang;}
		void setLabel(int label){Label=label;}
		void setNeigh(vector<int> neigh){Neigh=neigh;}
		void setVx(float vx){Vx=vx;}
		void setVy(float vy){Vy=vy;}
		
		//Calculate and run some internal things. --> MAKE THESE 4 METHODS OUTSIDE OF CLASS OBJECT
		void init();
		void updatePos(void);
		void storePos(void);
		void spreadSizes();

		//Return the variables for use.
		float getX(void){return X;}
		float getOldX(void){return OldX;}
		float getY(void){return Y;}
		float getOldY(void){return OldY;}
		float getVx(void){return Vx;}
		float getVy(void){return Vy;}
		float getAng(void){return Ang;}
		float getOldAng(void){return OldAng;}
		float getSize(void){return Size;}
		vector<int> getNeigh(void){return Neigh;}
		int getLabel(void){return Label;}


		particle();
};


particle::particle(void){
	X=0;
	Y=0;
	Ang=0;
	OldAng=0;
	Size=0;
	Label=0;
}
void particle::init(){ //For initialization.
	float dO=0;	//dO stores the distance from the particle to the origin.
	do {		//Assign a random X and Y to the particle within a box of size 2*L centered on the origin. 
		X=2*L*(float)rand()/(float)RAND_MAX-L;
		Y=2*L*(float)rand()/(float)RAND_MAX-L;
		dO=sqrtf(X*X+Y*Y);	//Calculate distance to the origin.
	} while (dO>L);			//If it is outside of the circle regenerate a random point untill you get one inside of the circle with radius L.
	Ang=2*M_PI*(float)rand()/(float)RAND_MAX;		//Random starting direction between 0 and 2*pi.
	OldAng=Ang;
	float rmu=log(RA/sqrtf((RS*RA*RS*RA/((float)(RA*RA)))+1));	//Mu and Sigma for generating size of the particle from a log normal distribution.
	float rsi=sqrtf(log(((RS*RA*RS*RA)/((float)(RA*RA)))+1));
	float u=(float)rand()/(float)RAND_MAX;				//u and v are random variables for generating a normal distribution.
	float v=(float)rand()/(float)RAND_MAX;
	Sizeeq=exp(rmu+rsi*sqrtf(-2*log(u))*cosf(2*M_PI*v));		//Generate an equilibrium size for particle from a log normal distribution. 
	Size=RA;							//Set the initial size to the average of the system. 
	Label=-1;							//Label everything -1 initially for later assignment to clusters.
}
void particle::updatePos(void){ 		//Move the position by the velocity amount. 
	X=X+Vx;
	Y=Y+Vy;
}
void particle::storePos(void){			//Save the current position, for use in calculating relative velocities for videos.
	OldX=X;
	OldY=Y;
}
void particle::spreadSizes(){			//Change the size by a small amount towards the equilibrium size, so after 100 changes it will be the right size.
	Size+=(Sizeeq-RA)/(100.0);
}

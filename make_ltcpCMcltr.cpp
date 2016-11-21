/* 	Calculate cluster related properties:
	1. csize: cluster counts with sizes in a specific region (1, >1, >5, >10, >30)
	2. csave: cluster size average (Define a cluster when # of atoms is more than 9)
	3. csind: cluster individual: all cluster sizes at several steps
	4.   lce: local chemical environment, percentage of solute atoms surrounding the vacancy during a period of time 
	5.   sro: short range order, 1-(zij/z/cj)
	6.   lro: long range order, abs[(-1)^i * sigma_i] / N

	Code structure:
	1. read_t0
	2. read_history
		a. update_cid(loop)
		b. sum csize
			(1) output csize
			(2) output csave
			(3) output csind
		c. cal sro
		d. cal lro
		e. cal lce

*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <array>
#include <algorithm>
#include <sstream>
#include <cmath>

using namespace std;

void error(int nexit, string errinfo, int nnum=0, double num1=0, double num2=0){
	cout << "\nError: ";
	cout << errinfo << " ";
	switch(nnum){
		case 0:	 cout << endl;			      break;
		case 1:  cout << num1 << endl; 		      break;
		case 2:  cout << num1 << " " << num2 << endl; break;
		default: cout << "!!!ERROR FUNCTION MALFUNCTION!!! WRONG NNUM!!!" << endl;
	}
	exit(1); 
}

int pbc(int x_, int nx_){ // Periodic Boundary Condition
	if	(x_<-nx_ || x_>=2*nx_)
		error(1, "(pbc) input x is out of bound", 2, x_, nx_);

	if	(x_<0) 		return (x_ + nx_);
	else if	(x_<nx_)	return  x_;
	else			return (x_ - nx_);
}

// parameters //
#define MAX_TYPE 10
#define CONST_CLTR 20 // The definition of minimum number of ltcps for a cluster
const double vbra[3][3]= {{-0.5,  0.5,  0.5}, { 0.5, -0.5,  0.5}, { 0.5,  0.5, -0.5}};
const int n1nbr= 8;
const int v1nbr[8][3]= {{ 1,  0,  0}, { 0,  1,  0}, { 0,  0,  1},
			{-1,  0,  0}, { 0, -1,  0}, { 0,  0, -1},
			{ 1,  1,  1}, {-1, -1, -1}	       };

const int nx= 128;
const int ny= 128;
const int nz= 128;

const int Ttype= 5; // targeted type: solute atom
// parameters //

// global variables
int def_cltr= CONST_CLTR;
double realtime= 0; // real time in the simulation (unit: s) 
long long int timestep= 0;

array <int, nx*ny*nz> states; // states[nx][ny][nz];

int cid= 0;
array <int, nx*ny*nz> id_cltr= {0}; // indicate the cluster id that a Ttype ltc(index) is labeled with (0 means not Ttype)

char name_in_ltcp[50];
char name_OUTPUT[50];
// global variables

//######################### functions #########################//
// Read files
void read_ltcp();
// Calculations
void update_cid(int x, int y, int z, vector <bool> & is_updated); // input one ltc point and renew cluster id around it
int  sum_csize(bool isdef= false);
//######################### functions #########################//

int main(int nArg, char *Arg[]){
	int t0cpu= time(0);

	// INPUT FILE NAMES
	if(nArg != 3){
		cout << "Make a ltcp file that the points are CM of clusters" << endl;
		cout << "Error: Parameter required: <name_in_ltcp> <name_OUTPUT>" << endl;
		exit(1);
	}
	strcpy(name_in_ltcp, Arg[1]);
	strcpy(name_OUTPUT, Arg[2]);
	// INPUT FILE NAMES
	
	cout << "Calculation is starting..." << endl;
	cout << "The TARGETED type is " << Ttype << endl;
	cout << "The system size is " << nx << " x " << ny << " x " << nz << endl;

	// CALCULATIONS
	read_ltcp();
	// CALCULATIONS

	int tfcpu= time(0);
	cout << "\nCalculation is completed!! Total CPU time: " << tfcpu - t0cpu << " secs" << endl;
}

void read_ltcp(){ // Reading t0.ltcp ////////////////////
    ifstream in_ltcp(name_in_ltcp);
	if(! in_ltcp.is_open()) error(1, "in_ltcp was not open"); // check

	int ntotal; in_ltcp >> ntotal; in_ltcp.ignore(); 
	if(ntotal != nx*ny*nz) error(1, "ntotal != nx*ny*nz", 2, ntotal, nx*ny*nz); // check
	
	string line2; getline(in_ltcp, line2);

    vector <int> states_temp(nx*ny*nz);
	for(int a=0; a<ntotal; a ++){
		int state_in, i, j, k, ix, iy, iz, srf, dump1, dump2;
		in_ltcp >> state_in >> i >> j >> k;
		if(state_in == 1 || state_in == -1) in_ltcp >> srf;
        else if(state_in == 2)              in_ltcp >> ix >> iy >> iz >> dump1 >> dump2;
        else                                in_ltcp >> ix >> iy >> iz; 
	
		if(in_ltcp.eof()) error(1, "reach end of file before finish reading all data");
		if(state_in > MAX_TYPE) error(1, "(in_t0) state input is larger than MAX_TYPE", 1, state_in); // check
		if(a != i*ny*nz+j*nz+k) error(1, "in_t0: ltcp inconsistent", 2, a, i*ny*nz+j*nz+k); // check
		states_temp[a]= state_in;
	}
	in_ltcp.close();


	// identify clusters and give them ids at t0
	for(int a=0; a<ntotal; a ++) states[a]= states_temp[a];

    vector <bool> is_updated(nx*ny*nz);
    is_updated= {false};
	for(int a=0; a<ntotal; a ++){
		int x= (int) (a/nz)/ny;
		int y= (int) (a/nz)%ny;
		int z= (int)  a%nz;

		update_cid(x, y, z, is_updated); // cltr
	}

	sum_csize(); // cltr
}
    
void update_cid(int x, int y, int z, vector <bool> & is_updated){ 
	int i_ltcp= x*ny*nz + y*nz + z;
	vector <int> i_ing;  //	the indexs which are being counted

	if(is_updated[i_ltcp]) goto skip_update;

	if(states[i_ltcp]==Ttype){
		cid ++;	
		
		i_ing.push_back(i_ltcp);
		id_cltr[i_ltcp]= cid;
		is_updated[i_ltcp]= true;

		for(int a=0; a<i_ing.size(); a++){ // !! be careful here: i_ing.size() increases during the iteration and automatically changes the condition
			int x1= (int) (i_ing.at(a)/nz)/ny;
			int y1= (int) (i_ing.at(a)/nz)%ny;
			int z1= (int)  i_ing.at(a)%nz;

			for(int b=0; b<n1nbr; b ++){
				int x2= pbc(x1+v1nbr[b][0], nx);
				int y2= pbc(y1+v1nbr[b][1], ny);
				int z2= pbc(z1+v1nbr[b][2], nz);
				int index= x2*ny*nz + y2*nz + z2;
	
				if(states[index]==Ttype){
					vector<int>::iterator it= find(i_ing.begin(), i_ing.end(), index);
					if(it==i_ing.end()){
						i_ing.push_back(index);
						id_cltr[index]= cid;
						is_updated[index]= true;
					}
				}
			}
		}
	}
	else id_cltr[i_ltcp]= 0;

skip_update:;
}

int sum_csize(bool isdef){
	int Ncsize[1000]= {0}; // # of clusters for single, l>1, l>5, l>10, l>30
	int Ncltr= 0;       // # of clusters, not atoms in clusters
	int sumNincltr= 0;  // total number of Ttype in clusters 
	double sumRadius= 0;// sum of radius of clusters
#define N_UCELL 2.0
#define PI 3.14159265359

	vector <int> N_in_cltr(cid+1); // for a certain cluster id, the number of counts
	vector <int> sum_x(cid+1);
	vector <int> sum_y(cid+1);
	vector <int> sum_z(cid+1);
    vector <int> firstx(cid+1); // the x position for first atom in the cltr; use for pbc check
    vector <int> firsty(cid+1);
    vector <int> firstz(cid+1);

    ofstream outCM(name_OUTPUT);

    for(int a=0; a<cid+1; a++){
        firstx[a]= -1;
        firsty[a]= -1;
        firstz[a]= -1;
    }

	for(int i=0; i<nx*ny*nz; i ++){
		int x= (int) (i/nz)/ny;
		int y= (int) (i/nz)%ny;
		int z= (int)  i%nz;
        int icltr= id_cltr[i];

		if(icltr != 0){
			N_in_cltr[icltr] ++;

            if(firstx[icltr]==-1) firstx[icltr]= x;
            if(firsty[icltr]==-1) firsty[icltr]= y;
            if(firstz[icltr]==-1) firstz[icltr]= z;

            if( (x-firstx[icltr]) >0.75*nx)     x -= nx;
            else if((x-firstx[icltr])<-0.75*nx) x += nx;
            if( (y-firsty[icltr]) >0.75*ny)     y -= ny;
            else if((y-firsty[icltr])<-0.75*ny) y += ny;
            if( (z-firstz[icltr]) >0.75*nz)     z -= nz;
            else if((z-firstz[icltr])<-0.75*nz) z += nz;

            sum_x[icltr] += x;
            sum_y[icltr] += y;
            sum_z[icltr] += z;
		}
	}

    int N= cid;
    for(int j=1; j<=cid; j ++) if(N_in_cltr[j]<def_cltr) N --;

    outCM << N << endl << endl;
    for(int j=1; j<=cid; j ++){
	    if(N_in_cltr[j]<def_cltr) continue;

        double cmx= sum_x[j]*1.0/N_in_cltr[j];
        double cmy= sum_y[j]*1.0/N_in_cltr[j];
        double cmz= sum_z[j]*1.0/N_in_cltr[j];

        if(abs(cmx/nx)>2) error(1, "cmx werid", 1, cmx);
        else if(cmx >=nx) cmx -= nx;
        else if(cmx <0)   cmx += nx;
        if(abs(cmy/ny)>2) error(1, "cmy werid", 1, cmy);
        else if(cmy >=ny) cmy -= ny;
        else if(cmy <0)   cmy += ny;
        if(abs(cmz/nz)>2) error(1, "cmz werid", 1, cmz);
        else if(cmz >=nz) cmz -= nz;
        else if(cmz <0)   cmz += nz;
    
        outCM << "6 " << cmx << " " << cmy << " " << cmz << endl;
    }
}

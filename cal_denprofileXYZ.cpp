#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <unordered_map>

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

// parameters //
const int nx=  64;
const int ny=  64;
const int nz=  64;

int Ttype1= -1;
int Ttype2=  3;
#define DX 4
#define AXIS x 
// parameters //

// global variables
const double vbra[3][3]= {{-0.5,  0.5,  0.5}, { 0.5, -0.5,  0.5}, { 0.5,  0.5, -0.5}};
long long int timestep; 
double realtime;
vector <int> states(nx*ny*nz);
vector <double> x_(nx*ny*nz);
vector <double> y_(nx*ny*nz);
vector <double> z_(nx*ny*nz);

int N_check= 0;

FILE * out_den; 

double lo=99999; 
double hi=-99999;
// global variables

//######################### functions #########################//
// Read files
void read_xyz(char ifname[]);
// Calculations
void cal_den();
//######################### functions #########################//

int main(int nArg, char *Arg[]){
	int t0cpu= time(0);

	cout << "Calculation is starting..." << endl;
	cout << "The system size is " << nx << " x " << ny << " x " << nz << endl;

	// READING FILES
	if(nArg != 2 && nArg != 3) error(0, "nArg must be 2 or 3\nUse den.exe <.xyz> <out>\n");
	read_xyz(Arg[1]);
	// READING FILES
	
	// OPEN OUTPUT FILES 
	char name[50]="out.", s_ts[50], s_time[50]; 
	if(nArg==3) strcpy(name, Arg[2]);
    else{
        sprintf(s_ts,   "%lld", timestep);
	    sprintf(s_time, "%.2e", realtime);
	    strcat(name, s_ts);
	    strcat(name, "_");
	    strcat(name, s_time);
    }
	
	out_den= fopen(name, "w");
	if(NULL==out_den) error(1, "out_den was not open");
	// OPEN OUTPUT FILES  
	
	// CALCULATIONS
	cal_den();
	// CALCULATIONS

	int tfcpu= time(0);
	cout << "\nCalculation is completed!! Total CPU time: " << tfcpu - t0cpu << " secs" << endl;
}

void read_xyz(char ifname[]){ // Reading .ltcp ////////////////////
	ifstream in_t0(ifname, ios::in);
	if(! in_t0.is_open()) error(1, "t0.xyz was not open"); // check

	int ntotal; in_t0 >> ntotal;
	if(ntotal != nx*ny*nz) error(1, "ntotal != nx*ny*nz", 2, ntotal, nx*ny*nz); // check

	char dump[50];
	in_t0 >> dump >> timestep >> realtime;
	cout << "the file is at " << timestep << " " << realtime << " (timestep time)" << endl;

	for(int a=0; a<ntotal; a ++){
        if(in_t0.eof()) error(1, "reach end of file before finish reading all data");

		int state_in, i, j, k, is_srf, ix, iy, iz, dir, head;
		in_t0 >> state_in >> i >> j >> k;

		if(state_in == 1 || state_in == -1) in_t0 >> is_srf;
        else if(state_in != 4)              in_t0 >> ix >> iy >> iz;

		if(state_in == 2 || state_in == -2 || state_in == 3) in_t0 >> dir >> head;
				
        double x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0];
		double y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1];
		double z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2];

		states[a]= state_in;
        x_[a]= x;
        y_[a]= y;
        z_[a]= z;
        if(AXIS>hi) hi= AXIS;
        if(AXIS<lo) lo= AXIS;
		if(Ttype1==state_in || Ttype2==state_in) N_check ++;
	}
	in_t0.close();
	cout << "t0.xyz file reading completed" << endl;
    cout << "hi, low: " << hi << " " << lo << endl;
}

void cal_den(){
	int N_Ttype= 0;

    unordered_map <double, int> NaAxis;
    unordered_map <double, int> NtAxis;
    for(int a=0; a<nx*ny*nz; a++){
        double x= x_[a];
        double y= y_[a];
        double z= z_[a];
        unordered_map <double, int>::iterator got = NaAxis.find(AXIS);
        if(got == NaAxis.end()){ NtAxis[AXIS]= 0; NaAxis[AXIS]= 0;}
        
        NaAxis.at(AXIS) ++;
	    if(Ttype1==states[a] || Ttype2==states[a]){
            N_Ttype ++;
            NtAxis.at(AXIS) ++;
        }
    }
	
    if(N_check != N_Ttype) error(2, "N Ttype inconsistent", N_check, N_Ttype);
	cout << "number of targeted type: " << N_Ttype << endl;
    
    for(int i=-nx*10; i<nx*10; i++){
        unordered_map <double, int>::iterator got = NaAxis.find(i);
        if(got != NaAxis.end()) fprintf(out_den, "%f %f\n", i/2.0, NtAxis[i]*1.0/NaAxis[i]);
    }
}

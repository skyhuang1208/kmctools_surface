#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
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

// parameters //
const int nx= 660;
const int ny= 270;
const int nz=  16;
const int nMlayer= 20;
// parameters //

// global variables
long long int timestep; 
double realtime;

FILE * out_den; 
// global variables

// functions //
void read_xyz(char ifname[]);
void cal_deg(double avg_conc, int den[nx]);
// functions //

int main(int nArg, char *Arg[]){
	// READING FILES
	if(nArg != 2 && nArg != 3) error(0, "nArg must be 2 or 3\nUse den.exe <.xyz> [file_name](optional)\n");
	
    // OPEN OUTPUT FILES 
	char name[50]="out.";
    if(nArg==3)	strcat(name, Arg[2]);
    else        strcat(name, "segdeg");
	
	out_den= fopen(name, "w");
	if(NULL==out_den) error(1, "out_den was not open");
	// OPEN OUTPUT FILES  
	
    read_xyz(Arg[1]);
	
    cout << "\nCalculation is completed!!" << endl;
}

void read_xyz(char ifname[]){ // Reading .ltcp ////////////////////
	ifstream in_t0(ifname, ios::in);
	if(! in_t0.is_open()) error(1, "t0.xyz was not open"); // check

	int ntotal;;
    while(in_t0 >> ntotal){
        int den[nx]= {0};

	    char dump[50];
	    in_t0 >> dump >> timestep >> realtime;
		if(strcmp("T:", dump) !=0) error(1, "(reading history) the format is incorrect"); // check

	    for(int a=0; a<ntotal; a ++){
            if(in_t0.eof()) error(1, "reach end of file before finish reading all data");

		    int ltcp; in_t0 >> ltcp;
		    int x= (int) (ltcp/nz)/ny;
            if(x<0 || x>=nx) error(1, "a out-of-bound x: (x nx)", 2, x, nx);

            den[x] ++;
	    }

        double avg_conc= double(ntotal)/(nx*ny*nz);
        cal_deg(avg_conc, den);
        cout << "T: " << timestep << " " << realtime << " (timestep time), avg conc= " << avg_conc << endl;
    }
	
    in_t0.close();
}

void cal_deg(double avg_conc, int den[nx]){
    int N_times= 0;
    double sum_deg= 0;
    
	for(int i=nMlayer; i<nx/2; i ++){
        double conc= double(den[i])/(ny*nz);
        if(conc < avg_conc) break;

        N_times ++;
        sum_deg += (conc - avg_conc);
    }
	for(int i=nx-nMlayer-1; i>nx/2; i --){
        double conc= double(den[i])/(ny*nz);
        if(conc < avg_conc) break;

        N_times ++;
        sum_deg += (conc - avg_conc);
    }

	double dis= 0.5*sqrt(3);
    cout << realtime << " " << sum_deg << endl;
	fprintf(out_den, "%f %f\n", realtime, dis*sum_deg);
}

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
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
const double vbra[3][3]= {{-0.5,  0.5,  0.5}, { 0.5, -0.5,  0.5}, { 0.5,  0.5, -0.5}};

double dr= 0.6;
double r_trunc= 128.0;

int type_A= -1;

const int nx= 64;
const int ny= 64;
const int nz= 64;
const double cx= 32.1761;
const double cy= 24.173;
const double cz= -4.50566;
// parameters //

void read_t0(char name[50], char name_out[50]);
double cal_dis(int dx, int dy, int dz);

int main(int nArg, char *Arg[]){
	int t0cpu= time(0);

	cout << "Calculation is starting..." << endl;
	cout << "The system size is " << nx << " x " << ny << " x " << nz << endl;

	// READ IMPORTED PARAMETERS
	if(nArg != 3){
		cout << "Calculate radial solute density profile" << endl;
		cout << "Error: Parameter required: <name_in_ltcp> <name_out>" << endl;
		exit(1);
	}

	// CALCULATIONS
	read_t0(Arg[1], Arg[2]);

	int tfcpu= time(0);
	cout << "\nCalculation is completed!! Total CPU time: " << tfcpu - t0cpu << " secs" << endl;
}

double cal_dis(int x, int y, int z){ // index pbc distance from the point (0, 0, 0)
	double dis= 9999999;
    for(int a=-2; a<=2; a++){
		for(int b=-2; b<=2; b++){
			for(int c=-2; c<=2; c++){
				int i= x + a*nx;
				int j= y + b*ny;
				int k= z + c*nz;
				double dx= cx - (i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0]);
				double dy= cy - (i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1]);
				double dz= cz - (i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2]);

				double d= sqrt(dx*dx + dy*dy + dz*dz);
                if(d<dis) dis= d;
	}}}

    return dis;
}

void read_t0(char name[50], char name_out[50]){ // Reading t0.ltcp
	cout << "\nNow reading t0.ltcp ...";
	
	ifstream in_t0(name, ios::in);
	if(! in_t0.is_open()) error(1, "in_t0 was not open"); // check

    int nall;
	in_t0 >> nall; in_t0.ignore();
	if(nall != nx*ny*nz) error(1, "ntotal != nx*ny*nz", 2, nall, nx*ny*nz); // check
	
	string line2; getline(in_t0, line2);

    vector <int> na( (int) (r_trunc/dr), 0);
    vector <int> nt( (int) (r_trunc/dr), 0);
    for(int a=0; a<nall; a ++){
		if(in_t0.eof()) error(1, "reach end of file before finish reading all data");
		
		int state_in;
		double i, j, k;
        int dump1, dump2, dump3, dump4, dump5;
		in_t0 >> state_in >> i >> j >> k;
        if(1==state_in || -1==state_in) in_t0 >> dump1;
        if(0==state_in ||  5==state_in) in_t0 >> dump1 >> dump2 >> dump3;
        if(2==state_in ||  3==state_in || -2==state_in) in_t0 >> dump1 >> dump2 >> dump3 >> dump4 >> dump5;

        if(i>(nx-1)) error(1, "point out of bound x", 1, i);
        if(j>(ny-1)) error(1, "point out of bound y", 1, j);
        if(k>(nz-1)) error(1, "point out of bound z", 1, k);

        double dis= cal_dis(i, j, k);

        na[(int) (dis/dr)] ++;
		if(type_A==state_in) nt[(int) (dis/dr)] ++;
	}
	cout << "t0.ltcp file reading completed" << endl;
	in_t0.close();
	
    FILE * OUT= fopen(name_out, "w");

    int ncheck= 0;
    for(int a=0; a<na.size(); a++){
        ncheck += na[a];
        if(na[a] != 0) fprintf(OUT, "%f %f %d %d\n", (a+0.5)*dr, nt[a]*1.0/na[a], nt[a], na[a]);
    }
    if(ncheck != nx*ny*nz) error(2, "number inconsistent", 2, ncheck, nx*ny*nz);
}

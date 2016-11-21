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
#define PI 3.14159265358979323846
const double vbra[3][3]= {{-0.5,  0.5,  0.5}, { 0.5, -0.5,  0.5}, { 0.5,  0.5, -0.5}};
const double sigma= 0.05; // standard deviation of the gaussian function: vibration distance 

double r0= 0.1;
double dr= 0.5;
double r_trunc= 64;
int nr; 

int type_A= 6;

const int nx= 128;
const int ny= 128;
const int nz= 128;
// parameters //

// global variables
int nall;
vector < vector<double> > Alist(3);

FILE * out;
char name_in_t0[20];
// global variables

//######################### functions #########################//
double cal_rho();
void read_t0();
void cal_rdf(double rho);
// Sub-functions for making cal lists
double cal_dis(int dx, int dy, int dz);
double gaussian(double x_);
// Sub-function for cal rdf
int pbc(int x_, int nx_);
//######################### functions #########################//

int main(int nArg, char *Arg[]){
	int t0cpu= time(0);

	cout << "Calculation is starting..." << endl;
	cout << "The system size is " << nx << " x " << ny << " x " << nz << endl;

	// READ IMPORTED PARAMETERS
	if(nArg != 3){
		cout << "Calculate RDF by using linked-cell method" << endl;
		cout << "WARNING!! have pre-assumed ltc are fixed BCC" << endl;
		cout << "Error: Parameter required: <name_in_ltcp> <name_out>" << endl;
		exit(1);
	}
	strcpy(name_in_t0, Arg[1]);
	cout << "\nImported parameters:" << endl << "dr= " << dr << ", truncation radius= " << r_trunc << endl;
	cout << "Atom type is: " << type_A << endl;

	// CAL PARAMETERS
	nr= r_trunc/dr + 10; // much larger than actually need 
	cout << "\nCalculated parameters: " << endl << "nr= " << nr << endl;

	// OPEN OUTPUT FILES  
	char name_nor[20];
	strcpy(name_nor, Arg[2]);
	out= fopen(name_nor, "w");
	if(NULL==out) error(1, "nor out file was not open");
	
	// CALCULATIONS
	read_t0();
	double rho= cal_rho(); cout << "Density is: " << rho << endl;
	cal_rdf(rho);

	int tfcpu= time(0);
	cout << "\nCalculation is completed!! Total CPU time: " << tfcpu - t0cpu << " secs" << endl;
}

double cal_rho(){ // calculate the density of the whole system
	double v1[3]= {nx*vbra[0][0], nx*vbra[0][1], nx*vbra[0][2]};
	double v2[3]= {ny*vbra[1][0], ny*vbra[1][1], ny*vbra[1][2]};
	double v3[3]= {nz*vbra[2][0], nz*vbra[2][1], nz*vbra[2][2]};

	double cross[3]= {v2[1]*v3[2]-v2[2]*v3[1], v2[2]*v3[0]-v2[0]*v3[2], v2[0]*v3[1]-v2[1]*v3[0]};

	if(r_trunc > nx*0.5*sqrt(3)) error(1, "truncation radius is larger than the system size");

	return nall / (v1[0]*cross[0] + v1[1]*cross[1] + v1[2]*cross[2]);
}

double cal_dis(int dx, int dy, int dz){ // index pbc distance from the point (0, 0, 0)
	double dis= 9999999;
    for(int a=-1; a<=1; a++){
		for(int b=-1; b<=1; b++){
			for(int c=-1; c<=1; c++){
				int da= dx + a*nx;
				int db= dy + b*ny;
				int dc= dz + c*nz;
				double dx2= da*vbra[0][0] + db*vbra[1][0] + dc*vbra[2][0];
				double dy2= da*vbra[0][1] + db*vbra[1][1] + dc*vbra[2][1];
				double dz2= da*vbra[0][2] + db*vbra[1][2] + dc*vbra[2][2];

				double d= sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
                if(d<dis) dis= d;
	}}}

    return dis;
}

double gaussian(double x_){
	if(abs(x_) > 3*sigma)	return 0;
	else                    return dr * ( exp(-x_*x_/2/sigma/sigma) / sqrt(2*PI*sigma*sigma) );
}

void read_t0(){ // Reading t0.ltcp
	cout << "\nNow reading t0.ltcp ...";
	
	ifstream in_t0(name_in_t0, ios::in);
	if(! in_t0.is_open()) error(1, "in_t0 was not open"); // check

	in_t0 >> nall; in_t0.ignore();
//	if(ntotal != nx*ny*nz) error(1, "ntotal != nx*ny*nz", 2, ntotal, nx*ny*nz); // check
	
	string line2; getline(in_t0, line2);

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

		if(type_A==state_in){
			Alist[0].push_back(i);
			Alist[1].push_back(j);
			Alist[2].push_back(k);
		}
	}
	cout << "t0.ltcp file reading completed" << endl;
	in_t0.close();
}

void cal_rdf(double rho){
	cout << "\nNow calculating rdf ..." << endl;

	double sum[nr]; // sum of all ltcp included within a certain shell
	for(int a=0; a < nr; a++) sum[a]= 0;

	for(int i=0; i<Alist[0].size(); i ++){
        double a= Alist[0].at(i);
        double b= Alist[1].at(i);
        double c= Alist[2].at(i);

	    for(int j=0; j<Alist[0].size(); j ++){
            double x= Alist[0].at(j);
			double y= Alist[1].at(j);
			double z= Alist[2].at(j);

            double dis= cal_dis(x-a, y-b, z-c);
            if(dis>r0 && dis<r_trunc){
                int nshell= (int) (dis/dr); // center of the shell
                for(int b= -3*sigma/dr-1; b<= 3*sigma/dr+1; b++){
                    double r_shell= (nshell+b+0.5)*dr;
			        double gvalue= gaussian(dis-r_shell);
                    sum[nshell+b] += gvalue;
                }
            }
		}
		
		if(0==i%1000) cout << i << "/" << Alist[0].size() << endl;
	}

	cout << "\nRDF calculations completed!! Now output to the file..." << endl; 

	for(int a=0; a < nr; a++){
		double r_print= 0.5*dr + a*dr;
		double r1= dr* a;
		double r2= dr*(a+1);
		double volume = 4*PI/3*(pow(r2,3) - pow(r1,3));

        fprintf(out, "%f %f\n", r_print, sum[a] / Alist[0].size() / volume / rho);
	}
}

int pbc(int x_, int nx_){ // Periodic Boundary Condition
	if	(x_<-nx_ || x_>=2*nx_)
		error(1, "(pbc) input x is out of bound", 2, x_, nx_);

	if	(x_<0) 		return (x_ + nx_);
	else if	(x_<nx_)	return  x_;
	else			return (x_ - nx_);
}

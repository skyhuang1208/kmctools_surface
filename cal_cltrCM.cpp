#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <array>
#include <algorithm>
#include <sstream>

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
#define MAX_TYPE 5
#define CONST_CLTR 3 // The definition of minimum number of ltcps for a cluster
const double vbra[3][3]= {{-0.5,  0.5,  0.5}, { 0.5, -0.5,  0.5}, { 0.5,  0.5, -0.5}};
const int n1nbr= 8;
const int v1nbr[8][3]= {{ 1,  0,  0}, { 0,  1,  0}, { 0,  0,  1},
			{-1,  0,  0}, { 0, -1,  0}, { 0,  0, -1},
			{ 1,  1,  1}, {-1, -1, -1}	       };

const int nx=  64;
const int ny=  64;
const int nz=  64;

const int Ttype= -1;
// parameters //

array <int, nx*ny*nz> states; // states[nx][ny][nz];
array <int, nx*ny*nz> id_cltr; // states[nx][ny][nz];
int cid= 0;
int ntt= 0;

//######################### functions #########################//
// Read files
void read_ltcp(char name_in_ltcp[50]);
void update_cid(int x, int y, int z, vector <bool> & is_updated); // input one ltc point and renew cluster id around it
void sum_csize();
//######################### functions #########################//

int main(int nArg, char *Arg[]){
	int t0cpu= time(0);

	// INPUT FILE NAMES
	if(nArg != 2){
		cout << "Calculate all CLUSTERING properties" << endl;
		cout << "Error: Parameter required: <name_in_ltcp>" << endl;
		exit(1);
	}
    char name[50];
	strcpy(name, Arg[1]);
	// INPUT FILE NAMES
	
	cout << "Calculation is starting..." << endl;
	cout << "The TARGETED type is " << Ttype << endl;
	cout << "The system size is " << nx << " x " << ny << " x " << nz << endl;

	// CALCULATIONS
	read_ltcp(name);
	// CALCULATIONS

	int tfcpu= time(0);
	cout << "\nCalculation is completed!! Total CPU time: " << tfcpu - t0cpu << " secs" << endl;
}

void read_ltcp(char name_in_ltcp[50]){ // Reading t0.ltcp ////////////////////
    ifstream in_ltcp(name_in_ltcp);
	if(! in_ltcp.is_open()) error(1, "in_ltcp was not open"); // check

	int ntotal; in_ltcp >> ntotal; in_ltcp.ignore(); 
	if(ntotal != nx*ny*nz) error(1, "ntotal != nx*ny*nz", 2, ntotal, nx*ny*nz); // check
	
	string line2; getline(in_ltcp, line2);

	for(int a=0; a<ntotal; a ++){
		int state_in, i, j, k, srf, ix, iy, iz, dir, head;
		
        in_ltcp >> state_in >> i >> j >> k;
		if(state_in == 1 || state_in == -1) in_ltcp >> srf;
        else if(state_in==4)                {}
        else{
            in_ltcp >> ix >> iy >> iz; 
            if(state_in != 0) in_ltcp >> head >> dir;
        }
		
		if(in_ltcp.eof()) error(1, "reach end of file before finish reading all data");
		if(state_in > MAX_TYPE) error(1, "(in_t0) state input is larger than MAX_TYPE", 1, state_in); // check
//		if(a != i*ny*nz+j*nz+k) error(1, "in_t0: ltcp inconsistent", 2, a, i*ny*nz+j*nz+k); // check
		int index= i*ny*nz + j*nz + k;
        states[index]= state_in;

		if(state_in==Ttype) ntt ++;
	}
	in_ltcp.close();
	cout << "t0.ltcp file reading completed; ntt= " << ntt << endl;

	// identify clusters and give them ids at t0
    int vcheck= 0;
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

void sum_csize(){
	int Ncltr= 0;       // # of clusters, not atoms in clusters
	int sumNincltr= 0;  // total number of Ttype in clusters 
#define N_UCELL 2.0
#define PI 3.14159265359

	vector <int> N_in_cltr(cid+1); // for a certain cluster id, the number of counts

	int N_check1= 0;
	for(int i=0; i<nx*ny*nz; i ++){
		if(id_cltr[i] != 0){
			N_in_cltr[id_cltr[i]] ++;
			N_check1 ++;
		}
	}

    vector <int> pickedID;
    vector <int> pickedSIZE;
    for(int a=0; a<(cid+1); a++){
        if(N_in_cltr[a]>100){
            pickedID.push_back(a);
            pickedSIZE.push_back(N_in_cltr[a]);
        }
    }

    vector <double> cx(pickedID.size(), 0);
    vector <double> cy(pickedID.size(), 0);
    vector <double> cz(pickedID.size(), 0);
    vector <int> x0(pickedID.size(), -1);
    vector <int> y0(pickedID.size(), -1);
    vector <int> z0(pickedID.size(), -1);
    vector <bool> sep(pickedID.size(), false); // if the cltr separate by PBC
    for(int i=0; i<nx*ny*nz; i++){
        for(int j= 0; j<pickedID.size(); j++){
            if(id_cltr[i]==pickedID[j]){
                int a= (int) (i/nz)/ny; // vcc position
	            int b= (int) (i/nz)%ny;
	            int c= (int)  i%nz;
	            
                if(-1==x0[j]){ x0[j]= a; y0[j]= b; z0[j]= c;}

                if(abs(a-x0[j])>nx/2 || abs(b-y0[j])>ny/2 || abs(c-z0[j])>nz/2) sep[j]= true;

                double x= a*vbra[0][0] + b*vbra[1][0] + c*vbra[2][0];
	            double y= a*vbra[0][1] + b*vbra[1][1] + c*vbra[2][1];
	            double z= a*vbra[0][2] + b*vbra[1][2] + c*vbra[2][2];

                cx[j] += x;
                cy[j] += y;
                cz[j] += z;
            }
        }
    }

    for(int a=0; a<pickedID.size(); a++){
        cx[a] /= pickedSIZE[a];
        cy[a] /= pickedSIZE[a];
        cz[a] /= pickedSIZE[a];

        cout << "cid, size, cx, cy, cz: " << pickedID[a] << " " << pickedSIZE[a] << " " << cx[a] << " " << cy[a] << " " << cz[a];
        if(sep[a]) cout << "\t<PBC separated!!>";
        cout << endl;
    }
}

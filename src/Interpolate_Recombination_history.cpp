//===========================================================================================================
// simple code to interpolate a given Xe-history as a function of redshift
//
// (JC, Feb 2016)
//===========================================================================================================

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include "Interpolate_Recombination_history.h"

using namespace std;

string Xefilename="./Xe-files/Xe_CosmoRec.dat";

bool file_is_loaded=0;
int glob_istart=0;

vector<double> zarr, Xearr;

//===========================================================================================================
double max(double a, double b){ return a>b ? a : b; }
double min(double a, double b){ return a<b ? a : b; }

//===========================================================================================================
// simple binomial search
//===========================================================================================================
void locate(const double xx[], unsigned long n, double x, unsigned long *j)
{
    // returns index j with xx[j]<= x; outside the boundaries it
    // will return the upper or lower index.
    unsigned long ju,jm,jl;
    int ascnd=(xx[n-1] > xx[0]);
    
    if(ascnd==1){
        if(x <= xx[0]){ *j=0; return; }
        if(x >= xx[n-1]){ *j=n-1; return; }
    }
    else{
        if(x >= xx[0]){ *j=0; return; }
        if(x <= xx[n-1]){ *j=n-1; return; }
    }
    
    jl=0;
    ju=n;
    while (ju-jl > 1) {
        jm=(ju+jl) >> 1;
        if ( (x >= xx[jm]) == ascnd)
            jl=jm;
        else
            ju=jm;
    }
    *j=jl;
    
    return;
}

//===========================================================================================================
// npol-1 is degree of the interpolating polynomial
//===========================================================================================================
void polint_routines(const double *xa, const double *ya,
                     int n, const double x,
                     double &y, double &dy)
{
    int i,m,ns=0;
    double den,dif,dift,ho,hp,w;
    
    vector<double> c(n),d(n);
    dif=fabs(x-xa[0]);
    for (i=0;i<n;i++) {
        if ((dift=fabs(x-xa[i])) < dif) {
            ns=i;
            dif=dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    y=ya[ns--];
    for (m=1;m<n;m++) {
        for (i=0;i<n-m;i++) {
            ho=xa[i]-x;
            hp=xa[i+m]-x;
            w=c[i+1]-d[i];
            if ((den=ho-hp) == 0.0)
            {
                cout << " polint error " << ho << " " << x << " " << hp << " "
                << i << " " << xa[i] << " " << ya[i] << " "
                << m << " " << xa[i+m] << " " << ya[i+m]
                << endl;
            }
            den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        }
        y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
    }
}

//===========================================================================================================
void polint(const double *xa, const double *ya, int na, const double x,
            int &istart, int npol, double *y, double *dy)
{
    // npol-1 is degree of the interpolating polynomial
    long unsigned int j=istart;
    locate(xa, na-1, x, &j);      // find index corresponding to x (start at i=istart)
    istart=j;
    
    if(istart==na-1){ *y=0.0; *dy=0.0; return; } // no point above x[j]!
    
    int npol_loc=(int)min(na, npol); // can only use as many as na points for interpolation
    int nlow=(npol_loc-1)/2; //, nup=npol_loc-1-nlow; // nlow+nup+1 == npol_loc!!!
    int ks=(int)min(max(istart-nlow, 0), na-1-npol_loc);  // defining region arround index
    
    polint_routines(&xa[ks], &ya[ks], npol_loc, x, *y, *dy);
    
    return;
}

//===========================================================================================================
void load_Xe_file()
{
    ifstream ifile(Xefilename.c_str());
    if(!ifile){ cerr << " Please check file: `" << Xefilename << "'" << endl; exit(1); }
    
    double dum;
    
    while(!ifile.eof())
    {
        ifile >> dum;
        zarr.push_back(log(dum));
        ifile >> dum;
        Xearr.push_back(log(dum));
        ifile >> dum;
    }
    return;
}

//===========================================================================================================
//
// z  == redshift
// Xe == Ne/NH where NH is the number density of Hydrogen nuclei NH ~ 1.9e-7 (1+z)^3 cm^-3 for Yp=0.24
//
//===========================================================================================================
void Compute_Xe(float *z, float *xe)
{
    if(!file_is_loaded)
    {
        load_Xe_file();
        file_is_loaded=1;
    }
    
    double lgz=log(min(*z, 8.0e+3));
    double y, dy;
    
    polint(&zarr[0], &Xearr[0], Xearr.size(), lgz, glob_istart, 3, &y, &dy);
    *xe = exp(y); 
}

//===========================================================================================================
//int main()
//{
//    cout << Compute_Xe(100.0) << " " << Compute_Xe(1100.0) << " " << Compute_Xe(1.0e+4) << endl;
//    return 0;
//}
//
//===========================================================================================================
//===========================================================================================================

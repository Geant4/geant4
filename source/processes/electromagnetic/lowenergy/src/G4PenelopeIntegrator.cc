//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
//
// Author:  Luciano Pandola (Luciano.Pandola@cern.ch)
//
// History:
// -----------
// 30 Nov 2002   LP        Created
//
// -------------------------------------------------------------------

#include "G4PenelopeIntegrator.hh"

// Constructor

template <class T,class F>
G4double G4PenelopeIntegrator<T,F>::Calculate(T& typeT, F f,
					      G4double LowPoint,G4double HighPoint,
					      G4double MaxError)
{
  const G4int npoints=10;
  const G4int ncallsmax=20000;
  const G4int nst=256;
  static G4double Abscissas[10] = {7.652651133497334e-02,2.2778585114164508e-01,3.7370608871541956e-01,
			    5.1086700195082710e-01,6.3605368072651503e-01,7.4633190646015079e-01,
			    8.3911697182221882e-01,9.1223442825132591e-01,9.6397192727791379e-01,
			    9.9312859918509492e-01};
  static G4double Weights[10] = {1.5275338713072585e-01,1.4917298647260375e-01,1.4209610931838205e-01,
			  1.3168863844917663e-01,1.1819453196151842e-01,1.0193011981724044e-01,
			  8.3276741576704749e-02,6.2672048334109064e-02,4.0601429800386941e-02,
			  1.7614007139152118e-02};

  //Error control
  G4double Ctol = G4std::min(G4std::max(MaxError,1e-13),1e-02);
  G4double Ptol = 0.01*Ctol;
  G4double Err=1e35;
  //Gauss integration from LowPoint to HighPoint
  G4double h,sumga,a,b,c,d;
  h=HighPoint-LowPoint;
  sumga=0.0;
  a=0.5*(HighPoint-LowPoint);
  b=0.5*(HighPoint+LowPoint);
  c=a*Abscissas[0];
  d=Weights[0]*((typeT.*f)(b+c)+(typeT.*f)(b-c));
  for (G4int i=2;i<=npoints;i++){
    c=a*Abscissas[i-1];
    d=d+Weights[i-1]*((typeT.*f)(b+c)+(typeT.*f)(b-c));
  }
  G4int icall = 2*npoints;
  G4int LH=1;
  G4double S[nst],x[nst],sn[nst],xrn[nst];
  S[0]=d*a;
  x[0]=LowPoint;

  //Adaptive bipartition scheme
  do{
   G4double h0=h;
   h=0.5*h; //bipartition
   G4double sumr=0;
   G4int LHN=0;
   G4double si,xa,xb,xc;
   for (i=1;i<=LH;i++){
    si=S[i-1];
    xa=x[i-1];
    xb=xa+h;
    xc=xa+h0;
    a=0.5*(xb-xa);
    b=0.5*(xb+xa);
    c=a*Abscissas[0];
    d=Weights[0]*((typeT.*f)(b+c)+(typeT.*f)(b-c));
    for (G4int j=1;j<npoints;j++){
      c=a*Abscissas[j];
      d=d+Weights[j]*((typeT.*f)(b+c)+(typeT.*f)(b-c));
    }    
    G4double s1,s2,s12;
    s1=d*a;
    a=0.5*(xc-xb);
    b=0.5*(xc+xb);
    c=a*Abscissas[0];
    d=Weights[0]*((typeT.*f)(b+c)+(typeT.*f)(b-c));
    for (j=1;j<npoints;j++){
      c=a*Abscissas[j];
      d=d+Weights[j]*((typeT.*f)(b+c)+(typeT.*f)(b-c));
    }    
    s2=d*a;
    icall=icall+4*npoints;
    s12=s1+s2;
    if (abs(s12-si)<G4std::max(Ptol*abs(s12),1e-35)){
      sumga=sumga+s12;
    }
    else
      {
	sumr=sumr+s12;
	LHN=LHN+2;
	sn[LHN-1]=s2;
	xrn[LHN-1]=xb;
	sn[LHN-2]=s1;
	xrn[LHN-2]=xa;
      }
    if (icall>ncallsmax || LHN>nst)
      {
	G4cout << "G4PenelopeIntegrator: " << G4endl;
	G4cout << "LowPoint: " << LowPoint << ", High Point: " << HighPoint << G4endl;
	G4cout << "Tolerance: " << MaxError << G4endl;
	G4cout << "Calls: " << icall << ", Integral: " << sumga << ", Error: " << Err << G4cout;
	G4cout << "Number of open subintervals: " << LHN << G4endl;
	G4cout << "WARNING: the required accuracy has not been attained" << G4endl;
	return sumga;
	  }
   }
   Err=abs(sumr)/G4std::max(abs(sumr+sumga),1e-35);
   if (Err < Ctol || LHN == 0) return sumga; //end of cycle
   LH=LHN;
   for (i=0;i<LH;i++){
     S[i]=sn[i];
     x[i]=xrn[i];
   }
  }while(Ctol < 1.0); //infinite loop
  return 0; //It should never get here
}
	
template <class T,class F>
G4double G4PenelopeIntegrator<T,F>::Calculate(T* ptrT, F f,
					      G4double LowPoint,G4double HighPoint,
					      G4double MaxError)
{
  return Calculate(*ptrT,f,LowPoint,HighPoint,MaxError);
} 
  
  

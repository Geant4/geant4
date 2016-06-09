//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
// Author:  Luciano Pandola (Luciano.Pandola@cern.ch)
//
// History:
// -----------
// 17 Feb 2003   LP        Created
// 17 Dec 2003   LP        Removed memory leak
// 17 Mar 2004   LP        Removed unnecessary calls to std::pow(a,b)
//
// -------------------------------------------------------------------

#include "G4PenelopeInterpolator.hh"
#include "G4DataVector.hh"

G4PenelopeInterpolator::G4PenelopeInterpolator(G4double* pX,G4double* pY,G4int nOfData,G4double S1,G4double SN) :
  a(0),b(0),c(0),d(0),x(0),y(0)
{
  // pX = X grid points (ascending order)
  // pY = corresponding function values
  // nOfData = number of points in the grid > 4
  // S1 and S2 = second derivatives at X[0] and X[nOfData-1]
  a = new G4DataVector;
  b = new G4DataVector;
  c = new G4DataVector;
  d = new G4DataVector;

  if (nOfData < 4 )
    {
      G4String excep = "Spline interpolation cannot be performed with less than 4 points";
      G4Exception(excep);
    }

  G4int N1=nOfData-1;
  G4int N2=nOfData-2;
  G4int i;
  G4int k=0;
  G4DataVector A(N1),B(N2),D(nOfData);//auxiliary arrays
  A.clear();
  B.clear();
  D.clear();
  for (i=0;i<N1;i++){
    if ((pX[i+1]-pX[i]) < 1.0e-13) 
      {
	G4String excep = "Spline x values not in increasing order";
	G4Exception(excep);
      }
    A.push_back(pX[i+1]-pX[i]);
    D.push_back((pY[i+1]-pY[i])/A[i]);
  }
 
  //Symmetric coefficient matrix
  for (i=0;i<N2;i++){
    B.push_back(2.0*(A[i]+A[i+1]));
    k=N1-i-1;
    D[k]=6.0*(D[k]-D[k-1]);
  }
  D[1]=D[1]-A[0]*S1;
  D[N1-1]=D[N1-1]-A[N1-1]*SN;
  
  //Gauss solution of the tridiagonal system
  for (i=1;i<N2;i++){
    G4double R=A[i]/B[i-1];
    B[i]=B[i]-R*A[i];
    D[i+1]=D[i+1]-R*D[i];
  }

  //The sigma coefficients are stored in array D
  D[N1-1]=D[N1-1]/B[N2-1];
  for (i=1;i<N2;i++){
    k=N1-i-1;
    D[k]=(D[k]-A[k]*D[k+1])/B[k-1];
  }
  D.push_back(SN);

  //Spline coefficients
  G4double SI1=S1;
  G4double SI=0,H=0,HI=0;
  G4double store=0;
  G4double help1=0;
  G4double help2=0;

  for (i=0;i<N1;i++){
    SI=SI1;
    SI1=D[i+1];
    H=A[i];
    HI=1.0/H;
    help1 = pX[i+1]*pX[i+1]*pX[i+1];
    help2 = pX[i]*pX[i]*pX[i];
    store=HI*(SI*help1-SI1*help2)/6.0+HI*(pY[i]*pX[i+1]-pY[i+1]*pX[i])+
      H*(SI1*pX[i]-SI*pX[i+1])/6.0;
    a->push_back(store);
    store=(HI/2.0)*(SI1*(pX[i]*pX[i])-SI*(pX[i+1]*pX[i+1]))+HI*(pY[i+1]-pY[i])+(H/6.0)*(SI-SI1);
    b->push_back(store);
    store=(HI/2.0)*(SI*pX[i+1]-SI1*pX[i]);
    c->push_back(store);
    store=(HI/6.0)*(SI1-SI);
    d->push_back(store);
  }
  //Natural cubic spline for x > x[nOfData-1]
  G4double FN=pY[nOfData-1];
  store=(*b)[N1-1]+pX[nOfData-1]*(2.0*(*c)[N1-1]+pX[nOfData-1]*3.0*(*d)[N1-1]);
  a->push_back(FN-pX[nOfData-1]*store);
  b->push_back(store);
  c->push_back(0.0);
  d->push_back(0.0);

  x = new G4DataVector;
  y = new G4DataVector;

  for (i=0;i<nOfData;i++){
    x->push_back(pX[i]);
    y->push_back(pY[i]);
  }
  return;
}

G4PenelopeInterpolator::~G4PenelopeInterpolator()
{
  delete a;
  delete b;
  delete c;
  delete d;
  delete x;
  delete y;
}

G4double G4PenelopeInterpolator::CubicSplineInterpolation(G4double xx)
{
  G4double interp=0;
  G4int index = FindBin(xx);
  interp=(*a)[index]+xx*((*b)[index]+xx*((*c)[index]+xx*(*d)[index]));
  return interp;
}

G4double G4PenelopeInterpolator::FirstDerivative(G4double xx)
{
  G4double interp=0;
  G4int index = FindBin(xx);
  interp=(*b)[index]+xx*((*c)[index]*2.0+xx*(*d)[index]*3.0);
  return interp;
}

G4double G4PenelopeInterpolator::CalculateMomentum(G4double UpperLimit,
							  G4int MomentumOrder)
{
  G4int i;
  G4int nOfData = (G4int) x->size();
  const G4double eps=1.0e-35;
  if (MomentumOrder < -1) G4Exception("Calculate Momentum: error 0");
  if (nOfData < 2) G4Exception("Calculate Momentum: error 1");
  if ((*x)[0]<0) G4Exception("Calculate Momentum: error 2");
  for (i=1;i<nOfData;i++)
    {
      if ((*x)[i]<0) G4Exception("Calculate Momentum: error 3");
      if ((*x)[i] < (*x)[i-1]) G4Exception ("Calculate Momentum: error 4");
    }

  G4double RMom=0.0;
  if (UpperLimit < (*x)[0]) return RMom;
  G4int iend=0;
  G4double xt=std::min(UpperLimit,(*x)[nOfData-1]);
  G4double x1,x2,y1,y2;
  G4double xtc,dx,dy,a1,b1,ds;
  for (i=0;i<(nOfData-1);i++){
    x1=std::max((*x)[i],eps);
    y1=(*y)[i];
    x2=std::max((*x)[i+1],eps);
    y2=(*y)[i+1];
    if (xt < x2) 
      {
      xtc=xt;
      iend=1;
      }
    else
      {
	xtc=x2;
      }
    dx=x2-x1;
    dy=y2-y1;
    if (std::abs(dx) > (1e-14*std::abs(dy))) 
      {
	b1=dy/dx;
	a1=y1-b1*x1;
	if (MomentumOrder == -1) 
	  {
	    ds=a1*std::log(xtc/x1)+b1*(xtc-x1);
	  }
	else
	  {
	    ds=a1*(std::pow(xtc,MomentumOrder+1)-std::pow(x1,MomentumOrder+1))/ ((G4double) (MomentumOrder+1))+
	      b1*(std::pow(xtc,MomentumOrder+2)-std::pow(x1,MomentumOrder+2))/((G4double) (MomentumOrder+2));
	  }
      }
    else
      {
	ds=0.5*(y1+y2)*std::pow((xtc-x1),MomentumOrder);
      }
    RMom += ds;
    if (iend != 0) return RMom;
  }
  return RMom;
}
		     
G4int G4PenelopeInterpolator::FindBin(G4double xx)
{
  //Finds the interval x[i],x[i+1] which contains the value xx

  G4int nbOfPoints=x->size();
  
  if (xx > (*x)[nbOfPoints-1])
    {
      return (nbOfPoints-1);
    }
  if (xx < (*x)[0])
    {
      return 0;
    }
  G4int i=0,i1=nbOfPoints-1;
  do{
    G4int it=(i+i1)/2;
    if (xx > (*x)[it]) 
      {
	i=it;
      }
    else
      {
	i1=it;
      }
  } while((i1-i) > 1);
  return i;
}



  

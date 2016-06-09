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
// $Id: G4GaussHermiteQ.cc,v 1.6 2005/03/15 19:11:35 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
#include "G4GaussHermiteQ.hh"


// ----------------------------------------------------------
//
// Constructor for Gauss-Hermite

G4GaussHermiteQ::G4GaussHermiteQ ( function pFunction, 
                                   G4int nHermite     ) 
   : G4VGaussianQuadrature(pFunction)
{
   const G4double tolerance = 1.0e-12 ;
   const G4int maxNumber = 12 ;
   
   G4int i=1, j=1, k=1 ;
   G4double newton=0.;
   G4double newton1=0.0, temp1=0.0, temp2=0.0, temp3=0.0, temp=0.0 ;
   G4double piInMinusQ = std::pow(pi,-0.25) ; // 1.0/std::sqrt(std::sqrt(pi)) ??

   fNumber = (nHermite +1)/2 ;
   fAbscissa = new G4double[fNumber] ;
   fWeight   = new G4double[fNumber] ;

   for(i=1;i<=fNumber;i++)
   {
      if(i == 1)
      {
         newton = std::sqrt((G4double)(2*nHermite + 1)) - 
                  1.85575001*std::pow((G4double)(2*nHermite + 1),-0.16666999) ;
      }
      else if(i == 2)
      {
         newton -= 1.14001*std::pow((G4double)nHermite,0.425999)/newton ;
      }
      else if(i == 3)
      {
         newton = 1.86002*newton - 0.86002*fAbscissa[0] ;
      }
      else if(i == 4)
      {
         newton = 1.91001*newton - 0.91001*fAbscissa[1] ;
      }
      else 
      {
         newton = 2.0*newton - fAbscissa[i - 3] ;
      }
      for(k=1;k<=maxNumber;k++)
      {
         temp1 = piInMinusQ ;
         temp2 = 0.0 ;
         for(j=1;j<=nHermite;j++)
         {
            temp3 = temp2 ;
            temp2 = temp1 ;
            temp1 = newton*std::sqrt(2.0/j)*temp2
                  - std::sqrt(((G4double)(j - 1))/j)*temp3 ;
         }
         temp = std::sqrt((G4double)2*nHermite)*temp2 ;
         newton1 = newton ;
         newton = newton1 - temp1/temp ;
         if(std::fabs(newton - newton1) <= tolerance) 
         {
            break ;
         }
      }
      if(k > maxNumber)
      {
         G4Exception("G4GaussHermiteQ::G4GaussHermiteQ()",
                     "OutOfRange", FatalException,
                     "Too many iterations in Gauss-Hermite constructor.") ;
      }
      fAbscissa[i-1] =  newton ;
      fWeight[i-1] = 2.0/(temp*temp) ;
   }
}


// ----------------------------------------------------------
//
// Gauss-Hermite method for integration of std::exp(-x*x)*nFunction(x)
// from minus infinity to plus infinity . 

G4double G4GaussHermiteQ::Integral() const 
{
   G4double integral = 0.0 ;
   for(G4int i=0;i<fNumber;i++)
   {
      integral += fWeight[i]*(fFunction(fAbscissa[i])
                + fFunction(-fAbscissa[i])) ;
   }
   return integral ;
}

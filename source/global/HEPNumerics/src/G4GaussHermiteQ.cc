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
// $Id: G4GaussHermiteQ.cc 67970 2013-03-13 10:10:06Z gcosmo $
//

#include "G4GaussHermiteQ.hh"
#include "G4PhysicalConstants.hh"

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
   G4double newton0=0.;
   G4double newton1=0.0, temp1=0.0, temp2=0.0, temp3=0.0, temp=0.0 ;
   G4double piInMinusQ = std::pow(pi,-0.25) ; // 1.0/std::sqrt(std::sqrt(pi)) ??

   fNumber = (nHermite +1)/2 ;
   fAbscissa = new G4double[fNumber] ;
   fWeight   = new G4double[fNumber] ;

   for(i=1;i<=fNumber;i++)
   {
      if(i == 1)
      {
         newton0 = std::sqrt((G4double)(2*nHermite + 1)) - 
                   1.85575001*std::pow((G4double)(2*nHermite + 1),-0.16666999) ;
      }
      else if(i == 2)
      {
         newton0 -= 1.14001*std::pow((G4double)nHermite,0.425999)/newton0 ;
      }
      else if(i == 3)
      {
         newton0 = 1.86002*newton0 - 0.86002*fAbscissa[0] ;
      }
      else if(i == 4)
      {
         newton0 = 1.91001*newton0 - 0.91001*fAbscissa[1] ;
      }
      else 
      {
         newton0 = 2.0*newton0 - fAbscissa[i - 3] ;
      }
      for(k=1;k<=maxNumber;k++)
      {
         temp1 = piInMinusQ ;
         temp2 = 0.0 ;
         for(j=1;j<=nHermite;j++)
         {
            temp3 = temp2 ;
            temp2 = temp1 ;
            temp1 = newton0*std::sqrt(2.0/j)*temp2
                  - std::sqrt(((G4double)(j - 1))/j)*temp3 ;
         }
         temp = std::sqrt((G4double)2*nHermite)*temp2 ;
         newton1 = newton0 ;
         newton0 = newton1 - temp1/temp ;
         if(std::fabs(newton0 - newton1) <= tolerance) 
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
      fAbscissa[i-1] =  newton0 ;
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

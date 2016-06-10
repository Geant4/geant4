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
// $Id: G4GaussLaguerreQ.cc 67970 2013-03-13 10:10:06Z gcosmo $
//
#include "G4GaussLaguerreQ.hh"



// ------------------------------------------------------------
//
// Constructor for Gauss-Laguerre quadrature method: integral from zero to
// infinity of std::pow(x,alpha)*std::exp(-x)*f(x).
// The value of nLaguerre sets the accuracy.
// The constructor creates arrays fAbscissa[0,..,nLaguerre-1] and 
// fWeight[0,..,nLaguerre-1] . 
//

G4GaussLaguerreQ::G4GaussLaguerreQ( function pFunction,
                                    G4double alpha,
                                    G4int nLaguerre      ) 
   : G4VGaussianQuadrature(pFunction)
{
   const G4double tolerance = 1.0e-10 ;
   const G4int maxNumber = 12 ;
   G4int i=1, k=1 ;
   G4double newton0=0.0, newton1=0.0,
            temp1=0.0, temp2=0.0, temp3=0.0, temp=0.0, cofi=0.0 ;

   fNumber = nLaguerre ;
   fAbscissa = new G4double[fNumber] ;
   fWeight   = new G4double[fNumber] ;
      
   for(i=1;i<=fNumber;i++)      // Loop over the desired roots
   {
      if(i == 1)
      {
         newton0 = (1.0 + alpha)*(3.0 + 0.92*alpha)
                 / (1.0 + 2.4*fNumber + 1.8*alpha) ;
      }
      else if(i == 2)
      {
         newton0 += (15.0 + 6.25*alpha)/(1.0 + 0.9*alpha + 2.5*fNumber) ;
      }
      else
      {
         cofi = i - 2 ;
         newton0 += ((1.0+2.55*cofi)/(1.9*cofi)
                    + 1.26*cofi*alpha/(1.0+3.5*cofi))
                    * (newton0 - fAbscissa[i-3])/(1.0 + 0.3*alpha) ;
      }
      for(k=1;k<=maxNumber;k++)
      {
         temp1 = 1.0 ;
         temp2 = 0.0 ;
         for(G4int j=1;j<=fNumber;j++)
         {
            temp3 = temp2 ;
            temp2 = temp1 ;
            temp1 = ((2*j - 1 + alpha - newton0)*temp2
                     - (j - 1 + alpha)*temp3)/j ;
         }
         temp = (fNumber*temp1 - (fNumber +alpha)*temp2)/newton0 ;
         newton1 = newton0 ;
         newton0 = newton1 - temp1/temp ;
         if(std::fabs(newton0 - newton1) <= tolerance) 
         {
            break ;
         }
      }
      if(k > maxNumber)
      {
         G4Exception("G4GaussLaguerreQ::G4GaussLaguerreQ()",
                     "OutOfRange", FatalException,
                     "Too many iterations in Gauss-Laguerre constructor") ;
      }
         
      fAbscissa[i-1] = newton0 ;
      fWeight[i-1] = -std::exp(GammaLogarithm(alpha + fNumber)
                   - GammaLogarithm((G4double)fNumber))/(temp*fNumber*temp2) ;
   }
}

// -----------------------------------------------------------------
//
// Gauss-Laguerre method for integration of
// std::pow(x,alpha)*std::exp(-x)*pFunction(x)
// from zero up to infinity. pFunction is evaluated in fNumber points
// for which fAbscissa[i] and fWeight[i] arrays were created in
// G4VGaussianQuadrature(double,int) constructor

G4double 
G4GaussLaguerreQ::Integral() const 
{
   G4double integral = 0.0 ;
   for(G4int i=0;i<fNumber;i++)
   {
      integral += fWeight[i]*fFunction(fAbscissa[i]) ;
   }
   return integral ;
}

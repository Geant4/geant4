// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              IONtest5.cc
//
// Version:		
// Date:		
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		17191/03/NL/LvH
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 6 October 2003, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "G4Bessel.hh"
#include "globals.hh"

#include <iomanip>
///////////////////////////////////////////////////////////////////////////////
//
int main()
{
//
//
// Provide simple banner at start of output.
//
  G4cout <<"COMMENCING IONtest5 ...." <<G4endl;
  G4cout <<G4endl;
  G4cout <<G4endl;
//
//
// Intantiate Bessel function object.
//
  G4Bessel besselFunction;
//
//
// The loop through from xmin to xmax and display results of calculation for
// standard (Numerical Recipes) algorithm.
//
  G4cout <<G4endl;
  G4cout <<G4endl;
  G4cout <<"EVALUATING STANDARD BESSEL FUNCTIONS " <<G4endl;
  G4cout <<G4endl;
  G4cout <<"              x"
         <<"          I0(x)          I1(x)          K0(x)          K1(x)"
         <<G4endl;
  G4cout <<"---------------"
         <<"------------------------------------------------------------"
         <<G4endl;
  G4cout.precision(9);
  G4cout.setf(std::ios::showpoint);
  for (G4double x=-100.0; x<=0.0; x+=0.1)
  {
    G4cout.setf(std::ios::fixed);
    G4cout <<std::setw(15) <<x;
    G4cout.setf(std::ios::scientific);
    G4cout <<std::setw(15) <<besselFunction.I0(x)
           <<std::setw(15) <<besselFunction.I1(x)
           <<std::setw(15) <<besselFunction.K0(x)
           <<std::setw(15) <<besselFunction.K1(x)
           <<G4endl;
  }
  for (G4double x=0.01; x<0.5; x+=0.01)
  {
    G4cout.setf(std::ios::fixed);
    G4cout <<std::setw(15) <<x;
    G4cout.setf(std::ios::scientific);
    G4cout <<std::setw(15) <<besselFunction.I0(x)
           <<std::setw(15) <<besselFunction.I1(x)
           <<std::setw(15) <<besselFunction.K0(x)
           <<std::setw(15) <<besselFunction.K1(x)
           <<G4endl;
  }
  for (G4double x=0.5; x<100.1; x+=0.1)
  {
    G4cout.setf(std::ios::fixed);
    G4cout <<std::setw(15) <<x;
    G4cout.setf(std::ios::scientific);
    G4cout <<std::setw(15) <<besselFunction.I0(x)
           <<std::setw(15) <<besselFunction.I1(x)
           <<std::setw(15) <<besselFunction.K0(x)
           <<std::setw(15) <<besselFunction.K1(x)
           <<G4endl;
  }
  G4cout <<"---------------"
         <<"------------------------------------------------------------"
         <<G4endl;
//
//
// The loop through from xmin to xmax and display results of calculation for
// more precise algorithm.
//
  G4cout <<G4endl;
  G4cout <<G4endl;
  G4cout <<"EVALUATING PRECISE BESSEL FUNCTIONS " <<G4endl;
  G4cout <<G4endl;
  G4cout <<"              x"
         <<"          I0(x)          I1(x)          K0(x)          K1(x)"
         <<G4endl;
  G4cout <<"---------------"
         <<"------------------------------------------------------------"
         <<G4endl;
  G4cout.precision(9);
  G4cout.setf(std::ios::showpoint);
  for (G4double x=-100.0; x<=0.0; x+=0.1)
  {
    G4cout.setf(std::ios::fixed);
    G4cout <<std::setw(15) <<x;
    G4cout.setf(std::ios::scientific);
    G4cout <<std::setw(15) <<besselFunction.pI0(x)
           <<std::setw(15) <<besselFunction.pI1(x)
           <<std::setw(15) <<besselFunction.pK0(x)
           <<std::setw(15) <<besselFunction.pK1(x)
           <<G4endl;
  }
  for (G4double x=0.01; x<0.5; x+=0.01)
  {
    G4cout.setf(std::ios::fixed);
    G4cout <<std::setw(15) <<x;
    G4cout.setf(std::ios::scientific);
    G4cout <<std::setw(15) <<besselFunction.pI0(x)
           <<std::setw(15) <<besselFunction.pI1(x)
           <<std::setw(15) <<besselFunction.pK0(x)
           <<std::setw(15) <<besselFunction.pK1(x)
           <<G4endl;
  }
  for (G4double x=0.5; x<100.1; x+=0.1)
  {
    G4cout.setf(std::ios::fixed);
    G4cout <<std::setw(15) <<x;
    G4cout.setf(std::ios::scientific);
    G4cout <<std::setw(15) <<besselFunction.pI0(x)
           <<std::setw(15) <<besselFunction.pI1(x)
           <<std::setw(15) <<besselFunction.pK0(x)
           <<std::setw(15) <<besselFunction.pK1(x)
           <<G4endl;
  }
  G4cout <<"---------------"
         <<"------------------------------------------------------------"
         <<G4endl;

  G4cout <<"TEST IONtest5 COMPLETE" <<G4endl;
  G4cout <<G4endl;
  G4cout <<G4endl;
  
  return 0;
}

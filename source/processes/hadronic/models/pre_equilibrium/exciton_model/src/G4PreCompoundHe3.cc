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
// $Id: G4PreCompoundHe3.cc 90337 2015-05-26 08:34:27Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PreCompoundHe3
//
// Author:         V.Lara
//
// Modified:  
// 21.08.2008 J. M. Quesada add choice of options  
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup
// 05.07.2013 J.M. Quesada FactorialFactor fixed
//
 
#include "G4PreCompoundHe3.hh"
#include "G4SystemOfUnits.hh"
#include "G4He3.hh"

G4PreCompoundHe3::G4PreCompoundHe3()
  : G4PreCompoundIon(G4He3::He3(), &theHe3CoulombBarrier)
{}

G4PreCompoundHe3::~G4PreCompoundHe3()
{}

G4double G4PreCompoundHe3::FactorialFactor(G4int N, G4int P) const
{
  return G4double((N-3)*(P-2)*(N-2)*(P-1)*(N-1)*P)/12.0; 
}
  
G4double G4PreCompoundHe3::CoalescenceFactor(G4int A) const
{
  return 243.0/G4double(A*A);
}    

G4double G4PreCompoundHe3::GetRj(G4int nParticles, G4int nCharged) const
{
  G4double rj = 0.0;
  if(nCharged >=2 && (nParticles-nCharged) >= 1) {
    G4double denominator = G4double(nParticles*(nParticles-1)*(nParticles-2));
    rj = G4double(3*nCharged*(nCharged-1)*(nParticles-nCharged))/denominator;  
  }
  return rj;
}

G4double G4PreCompoundHe3::GetAlpha() const
{
  G4double C = 0.0;
  if (theFragZ <= 30) 
    {
      C = 0.10;
    }
  else if (theFragZ <= 50) 
    {
      C = 0.1 - (theFragZ - 30)*0.001;
    } 
  else if (theFragZ < 70) 
    {
      C = 0.08 - (theFragZ - 50)*0.001;
    }
  else 
    {
      C = 0.06;
    }
  return 1.0 + C*(4.0/3.0);
}


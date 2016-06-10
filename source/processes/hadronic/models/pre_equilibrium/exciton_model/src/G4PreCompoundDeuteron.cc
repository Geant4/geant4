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
// $Id: G4PreCompoundDeuteron.cc 90337 2015-05-26 08:34:27Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PreCompoundDeuteron
//
// Author:         V.Lara
//
// Modified:  
// 21.08.2008 J. M. Quesada add choice of options  
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup
//

#include "G4PreCompoundDeuteron.hh"
#include "G4SystemOfUnits.hh"
#include "G4Deuteron.hh"

G4PreCompoundDeuteron::G4PreCompoundDeuteron()
  : G4PreCompoundIon(G4Deuteron::Deuteron(), &theDeuteronCoulombBarrier)
{}

G4PreCompoundDeuteron::~G4PreCompoundDeuteron()
{}
 
G4double G4PreCompoundDeuteron::FactorialFactor(G4int N, G4int P) const
{
  return G4double((N-1)*(N-2)*(P-1)*P)*0.5;
}
  
G4double G4PreCompoundDeuteron::CoalescenceFactor(G4int A) const
{
  return 16.0/G4double(A);
}    

G4double G4PreCompoundDeuteron::GetRj(G4int nParticles, G4int nCharged) const
{
  G4double rj = 0.0;
  if(nCharged >=1 && (nParticles-nCharged) >=1) {
    G4double denominator = G4double(nParticles*(nParticles-1));
    rj = 2*nCharged*(nParticles-nCharged)/denominator; 
  }
  return rj;
}

G4double G4PreCompoundDeuteron::GetAlpha() const
{
  G4double C = 0.0;
  if (theFragZ >= 70) 
    {
      C = 0.10;
    } 
  else 
    {
      C = ((((0.15417e-06*theFragZ) - 0.29875e-04)*theFragZ 
	    + 0.21071e-02)*theFragZ - 0.66612e-01)*theFragZ + 0.98375; 
    }
  return 1.0 + 0.5*C;
}

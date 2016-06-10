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
// $Id: G4HETCProton.cc 90337 2015-05-26 08:34:27Z gcosmo $
//
// by V. Lara
//
// Modified:
// 23.08.2010 V.Ivanchenko general cleanup, move constructor and destructor 
//            the source, use G4Pow

#include "G4HETCProton.hh"
#include "G4Proton.hh"

G4HETCProton::G4HETCProton() 
  : G4HETCChargedFragment(G4Proton::Proton(), &theProtonCoulombBarrier)
{}

G4HETCProton::~G4HETCProton() 
{}

G4double G4HETCProton::GetAlpha() const
{
  G4double C = 0.0;
  if (theResZ >= 70) 
    {
      C = 0.10;
    } 
  else 
    {
      C = ((((0.15417e-06*theResZ) - 0.29875e-04)*theResZ 
	    + 0.21071e-02)*theResZ - 0.66612e-01)*theResZ + 0.98375;
    }
  return 1.0 + C;
}
  
G4double G4HETCProton::GetBeta() const
{
  return -theCoulombBarrier;
}
  
G4double G4HETCProton::GetSpinFactor() const
{
  // 2s+1
  return 2.0;
}

G4double G4HETCProton::K(const G4Fragment & aFragment)
{
  // Number of protons in emitted fragment
  G4int Pa = theZ;

  G4double r = G4double(theResZ)/G4double(theResA);

  G4int P = aFragment.GetNumberOfParticles();
  G4int H = aFragment.GetNumberOfHoles();

  G4double result = 0.0;
  if (P > 0)
    {
      result = (H*r + Pa)/(P*r);
    }

  return std::max(0.0,result);
}

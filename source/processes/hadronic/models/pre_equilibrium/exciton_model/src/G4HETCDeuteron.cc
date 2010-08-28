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
// $Id: G4HETCDeuteron.cc,v 1.3 2010-08-28 15:16:55 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara
//
// Modified:
// 23.08.2010 V.Ivanchenko general cleanup, move constructor and destructor 
//            the source, use G4Pow

#include "G4HETCDeuteron.hh"
#include "G4Deuteron.hh"

G4HETCDeuteron::G4HETCDeuteron() 
  : G4HETCChargedFragment(G4Deuteron::Deuteron(), &theDeuteronCoulombBarrier)
{}

G4HETCDeuteron::~G4HETCDeuteron() 
{}

G4double G4HETCDeuteron::GetAlpha()
{
  G4double C = 0.0;
  G4int aZ = GetZ() + GetRestZ();
  if (aZ >= 70) 
    {
      C = 0.10;
    } 
  else 
    {
      C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375; 
    }
  return 1.0 + C/2.0;
}
  
G4double G4HETCDeuteron::GetBeta()
{
  return -GetCoulombBarrier();
}

G4double G4HETCDeuteron::GetSpinFactor()
{
  // 2s+1
  return 3.0;
}

G4double G4HETCDeuteron::K(const G4Fragment & aFragment)
{
  // Number of protons in emitted fragment
  G4int Pa = GetZ();
  // Number of neutrons in emitted fragment 
  G4int Na = GetA() - Pa;

  G4int TargetZ = GetRestZ();
  G4int TargetA = GetRestA();
  G4double r = G4double(TargetZ)/G4double(TargetA);
  
  G4int P = aFragment.GetNumberOfParticles();
  G4int H = aFragment.GetNumberOfHoles();

  G4double result = 0.0;
  if (P > 1)
    {
      result = 2.0* (H*(H-1.0)*r*(r-1.0)+H*(Na*r+Pa*(1.0-r)) + Pa*Na)/(P*(P-1.0));

      result /= r*(1.0 - r);
    }
  return std::max(0.0,result);
}

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
// $Id$
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#include "G4GEMCoulombBarrierHE.hh"
#include "G4HadronicException.hh"
#include "G4Pow.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4GEMCoulombBarrierHE::G4GEMCoulombBarrierHE(G4int anA, G4int aZ) :
  G4VCoulombBarrier(anA,aZ) 
{}

G4GEMCoulombBarrierHE::~G4GEMCoulombBarrierHE() 
{}

G4double G4GEMCoulombBarrierHE::GetCoulombBarrier(G4int ARes, G4int ZRes, G4double U) const 
  // Calculation of Coulomb potential energy (barrier) for outgoing fragment
{
  G4double Barrier = 0.0;
  if (ZRes > ARes || ARes < 1) {
    G4cout << "G4GEMCoulombBarrierHE::GetCoulombBarrier: "
	   << "Wrong values for "
	   << "residual nucleus A = " << ARes << " "
	   << "and residual nucleus Z = " << ZRes << G4endl;
    throw G4HadronicException(__FILE__, __LINE__,"FATAL error");
  }
  if (GetZ() == 0) {
    Barrier = 0.0;   // If there is no charge there is neither barrier
  } else {
    G4double CompoundRadius = CalcCompoundRadius(ARes);
    Barrier = ( elm_coupling * GetZ() * ZRes)/(CompoundRadius+3.75*fermi);
    
    // Barrier penetration coeficient
    //    G4double K = BarrierPenetrationFactor(ZRes);
    //    Barrier *= K;
    
    Barrier /= (1.0 + std::sqrt(U/static_cast<G4double>(2*ARes)));
  }
  return Barrier;
}


G4double G4GEMCoulombBarrierHE::CalcCompoundRadius(G4int ARes) const
{
  G4Pow* g4pow = G4Pow::GetInstance();
  G4double AresOneThird = g4pow->Z13(ARes);
  G4double AejectOneThird = g4pow->Z13(GetA());

  G4double Result = 1.12*(AresOneThird + AejectOneThird) - 
    0.86*(AresOneThird+AejectOneThird)/(AresOneThird*AejectOneThird);

  return Result*fermi;
}



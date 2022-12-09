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
// J. M. Quesada (July 2009):  New class based on G4GEMCoulombBarrierHE
// Coded strictly according to Furihata's GEM paper 
// NEW:effective decrease  of barrier with E* (Barashenkov) has been added
//
#include "G4GEMCoulombBarrier.hh"
#include "G4HadronicException.hh"
#include "G4Pow.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4GEMCoulombBarrier::G4GEMCoulombBarrier(G4int anA, G4int aZ) :
  G4CoulombBarrier(anA, aZ) 
{
  AejectOneThird = g4calc->Z13(anA);
}

G4double G4GEMCoulombBarrier::GetCoulombBarrier(G4int ARes, G4int ZRes, 
                                                G4double U) const 
{
  // Calculation of Coulomb potential energy (barrier) for outgoing fragment
  G4double Barrier = 0.0;
  if (theZ > 0 && ZRes > 0) {

    G4double CompoundRadius = CalcCompoundRadius(ARes);
    Barrier = CLHEP::elm_coupling * (theZ * ZRes)/CompoundRadius;
      
    // Barrier penetration coeficient
    if(theA <= 4) { Barrier *= BarrierPenetrationFactor(ZRes); }
  
    //JMQ 200709 effective decrease  of barrier with E* (Barashenkov)
    // (not inclued in original Furihata's formulation)
    Barrier /= (1.0 + std::sqrt(U/((2*ARes)*CLHEP::MeV)));
  }
  return Barrier;
}

G4double G4GEMCoulombBarrier::CalcCompoundRadius(G4int ARes) const
{      
  G4double AresOneThird = g4calc->Z13(ARes);

  G4double Result = 0.0;
  if(theA == 1){
    Result = 1.7* AresOneThird;

  } else if (theA <= 4){
    Result = 1.7* AresOneThird + 1.2;

  } else {
    Result = 1.12*(AresOneThird + AejectOneThird) - 
      0.86*(AresOneThird+AejectOneThird)/(AresOneThird*AejectOneThird)+3.75;
  }
  return Result*CLHEP::fermi;
}



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
// $Id: G4MuonicAtomHelper.cc 96797 2016-05-09 10:13:42Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	History: first implementation, 1 July 16 K.Lynch
// ---------------------------------------------------------------

#include "G4MuonicAtomHelper.hh"

G4MuonicAtom* 
G4MuonicAtomHelper::ConstructMuonicAtom(G4String name, G4int encoding, G4Ions const* baseion){

  // mass calculation is missing binding energy correction
  //auto const mass = G4MuonMinus::Definition()->GetPDGMass() + baseion->GetPDGMass();
  auto const mass = 0.1056583715*CLHEP::GeV + baseion->GetPDGMass();
  // what should static charge be?  for G4Ions, it is Z ... should it
  // be Z-1 here (since there will always be a muon attached), or Z?
  auto const charge = baseion->GetPDGCharge(); 
  
  auto muatom = new G4MuonicAtom(name, mass, 0.0, charge, 
				 baseion->GetPDGiSpin(),
				 baseion->GetPDGiParity(),
				 baseion->GetPDGiConjugation(),
				 baseion->GetPDGiIsospin(),
				 baseion->GetPDGiIsospin3(),
				 baseion->GetPDGiGParity(),
				 baseion->GetParticleType(),
				 baseion->GetLeptonNumber(),
				 baseion->GetBaryonNumber(),
				 encoding,
				 baseion->GetPDGStable(),
				 // TODO: this is wrong ... need to
				 // get correct muonic atom lifetime.
				 // but see next comment
				 baseion->GetPDGLifeTime(),
				 // TODO: this is _definitely_ wrong:
				 // we need to add muonic atom decay
				 // channels ... do we try to support
				 // other modes of decay, or do we
				 // just assume for now that we don't
				 // need to worry about radioactive
				 // decay and the like?
				 baseion->GetDecayTable(),
				 baseion->IsShortLived(),
				 baseion->GetParticleSubType());

  muatom->SetPDGMagneticMoment(baseion->GetPDGMagneticMoment());
  return muatom;
}


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
// $Id: G4FinalStateChargeDecrease.hh,v 1.2 2008-07-14 20:47:34 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4FINALSTATECHARGEDECREASE_HH
#define G4FINALSTATECHARGEDECREASE_HH 1
 
#include "G4FinalStateProduct.hh"
#include "G4CrossSectionChargeDecreasePartial.hh"

class G4FinalStateChargeDecrease
{
 public:
   
   G4FinalStateChargeDecrease();
   
   ~G4FinalStateChargeDecrease();
   
   const G4FinalStateProduct& GenerateFinalState(const G4Track& track, const G4Step& step);
   
 private:
   
   G4double lowEnergyLimit;
   G4double highEnergyLimit;
   G4FinalStateProduct product;
   
   G4CrossSectionChargeDecreasePartial cross;  
 
   G4int NumberOfFinalStates(G4ParticleDefinition* particleDefinition, G4int finalStateIndex);
   G4ParticleDefinition* OutgoingParticleDefinition(G4ParticleDefinition* particleDefinition, G4int finalStateIndex);
   G4double WaterBindingEnergyConstant(G4ParticleDefinition* particleDefinition, G4int finalStateIndex);
   G4double OutgoingParticleBindingEnergyConstant(G4ParticleDefinition* particleDefinition, G4int finalStateIndex);

};

#endif

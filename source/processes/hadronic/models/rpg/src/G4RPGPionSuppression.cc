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
// $Id: G4RPGPionSuppression.cc 79697 2014-03-12 13:10:09Z gcosmo $
//
 
#include <iostream>
#include <signal.h>

#include "G4RPGPionSuppression.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4HadReentrentException.hh"

G4RPGPionSuppression::G4RPGPionSuppression()
  : G4RPGReaction() {}


G4bool G4RPGPionSuppression::
ReactionStage(const G4HadProjectile* /*originalIncident*/,
              G4ReactionProduct& modifiedOriginal,
              G4bool& incidentHasChanged,
              const G4DynamicParticle* /*originalTarget*/,
              G4ReactionProduct& targetParticle,
              G4bool& targetHasChanged,
              const G4Nucleus& targetNucleus,
              G4ReactionProduct& currentParticle,
              G4FastVector<G4ReactionProduct,256>& vec,
              G4int& vecLen,
              G4bool /*leadFlag*/,
              G4ReactionProduct& /*leadingStrangeParticle*/)
{
  // This code was originally in the fortran code TWOCLU
  //
  // Suppress charged pions, for various reasons
  //
  G4double mOriginal = modifiedOriginal.GetMass()/GeV;
  G4double etOriginal = modifiedOriginal.GetTotalEnergy()/GeV;
  G4double targetMass = targetParticle.GetDefinition()->GetPDGMass()/GeV;
  G4double cmEnergy = std::sqrt( mOriginal*mOriginal + targetMass*targetMass +
		      	   2.0*targetMass*etOriginal ); 
  G4double eAvailable = cmEnergy - mOriginal - targetMass;
  for (G4int i = 0; i < vecLen; i++) eAvailable -= vec[i]->GetMass()/GeV;

  const G4double atomicWeight = targetNucleus.GetA_asInt();
  const G4double atomicNumber = targetNucleus.GetZ_asInt();
  const G4double pOriginal = modifiedOriginal.GetTotalMomentum()/GeV;
    
  G4ParticleDefinition *aPiMinus = G4PionMinus::PionMinus();
  G4ParticleDefinition *aPiPlus = G4PionPlus::PionPlus();
  G4ParticleDefinition* aPiZero = G4PionZero::PionZero();
  G4ParticleDefinition *aProton = G4Proton::Proton();
  G4ParticleDefinition *aNeutron = G4Neutron::Neutron();
  G4double piMass = aPiPlus->GetPDGMass()/GeV;
  G4double nucleonMass = aNeutron->GetPDGMass()/GeV;
    
  const G4bool antiTest =
    modifiedOriginal.GetDefinition() != G4AntiProton::AntiProton() &&
    modifiedOriginal.GetDefinition() != G4AntiNeutron::AntiNeutron() &&
    modifiedOriginal.GetDefinition() != G4AntiLambda::AntiLambda() &&
    modifiedOriginal.GetDefinition() != G4AntiSigmaPlus::AntiSigmaPlus() &&
    modifiedOriginal.GetDefinition() != G4AntiSigmaMinus::AntiSigmaMinus() &&
    modifiedOriginal.GetDefinition() != G4AntiXiZero::AntiXiZero() &&
    modifiedOriginal.GetDefinition() != G4AntiXiMinus::AntiXiMinus() &&
    modifiedOriginal.GetDefinition() != G4AntiOmegaMinus::AntiOmegaMinus();

  if( antiTest && (
        currentParticle.GetDefinition() == aPiPlus ||
        currentParticle.GetDefinition() == aPiZero ||
        currentParticle.GetDefinition() == aPiMinus ) &&
      ( G4UniformRand() <= (10.0-pOriginal)/6.0 ) &&
      ( G4UniformRand() <= atomicWeight/300.0 ) )
  {
    if (eAvailable > nucleonMass - piMass) {
      if( G4UniformRand() > atomicNumber/atomicWeight )
        currentParticle.SetDefinitionAndUpdateE( aNeutron );
      else
        currentParticle.SetDefinitionAndUpdateE( aProton );
      incidentHasChanged = true;
    }
  }

  if( antiTest && (
        targetParticle.GetDefinition() == aPiPlus ||
        targetParticle.GetDefinition() == aPiZero ||
        targetParticle.GetDefinition() == aPiMinus ) &&
      ( G4UniformRand() <= (10.0-pOriginal)/6.0 ) &&
      ( G4UniformRand() <= atomicWeight/300.0 ) )
  {
    if (eAvailable > nucleonMass - piMass) {
      if( G4UniformRand() > atomicNumber/atomicWeight )
        targetParticle.SetDefinitionAndUpdateE( aNeutron );
      else
        targetParticle.SetDefinitionAndUpdateE( aProton );
      targetHasChanged = true;
    }
  }

  for( G4int i=0; i<vecLen; ++i )
  {
    if( antiTest && (
          vec[i]->GetDefinition() == aPiPlus ||
          vec[i]->GetDefinition() == aPiZero ||
          vec[i]->GetDefinition() == aPiMinus ) &&
        ( G4UniformRand() <= (10.0-pOriginal)/6.0 ) &&
        ( G4UniformRand() <= atomicWeight/300.0 ) )
    {
      if (eAvailable > nucleonMass - piMass) {
        if( G4UniformRand() > atomicNumber/atomicWeight )
          vec[i]->SetDefinitionAndUpdateE( aNeutron );
        else
          vec[i]->SetDefinitionAndUpdateE( aProton );
      }
    }
  }

  return true;
}

 
 /* end of file */

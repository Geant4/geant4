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

#include "Tst15SteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4IsoParticleChange.hh"
#include "G4InelasticInteraction.hh"
#include "G4ios.hh"
#include <iomanip>

Tst15SteppingAction::Tst15SteppingAction()
{}

Tst15SteppingAction::~Tst15SteppingAction()
{}

void Tst15SteppingAction::UserSteppingAction(const G4Step* theStep)
{
  // Is track alive?
  G4Track* theTrack = theStep->GetTrack();
  if (theTrack->GetTrackStatus() != fAlive) {
    // No, print isotope production information
    G4String procName =
      theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    if (procName == "NeutronInelastic" || procName == "ProtonInelastic") {
      G4IsoParticleChange* isotopeProdInfo =
        G4InelasticInteraction::GetIsotopeProductionInfo();
      if (isotopeProdInfo) {
        G4cout << " Production model: " << isotopeProdInfo->GetProducer() << G4endl;
        G4cout << " target nucleus (A,Z) : "
               << isotopeProdInfo->GetMotherNucleus().GetA_asInt() << " , "
               << isotopeProdInfo->GetMotherNucleus().GetZ_asInt() << G4endl; 
        G4cout << " produced isotope: " << isotopeProdInfo->GetIsotope() << G4endl;
        G4cout << " production time: " << isotopeProdInfo->GetProductionTime() << G4endl;
        G4cout << " parent particle: "
               << isotopeProdInfo->GetParentParticle()->GetParticleName() << G4endl;
      } 
    }
  }
}


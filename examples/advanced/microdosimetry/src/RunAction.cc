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
// -------------------------------------------------------------------
// $Id$
// -------------------------------------------------------------------

#include "RunAction.hh"
#include "G4Run.hh"
#include "TrackingAction.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::RunAction(DetectorConstruction* det, HistoManager* his, TrackingAction* trackingAction)
    :Detector(det),Histo(his),TrackingAct(trackingAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run*)
{  
  // Histograms
  Histo->book();
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrintNParticles(std::map<const G4ParticleDefinition*, int>& container)
{
    std::map<const G4ParticleDefinition*, int>::iterator it;
    for(it = container.begin() ;
        it != container.end(); it ++)
    {
        G4cout << "N " << it->first->GetParticleName() << " : " << it->second << G4endl;
    }
}

void RunAction::EndOfRunAction(const G4Run*)
{
  //save histograms      
  Histo->save();

  std::map<const G4ParticleDefinition*, int>&  particlesCreatedInWorld = TrackingAct->GetNParticlesCreatedInWorld();
  G4cout << "Number and type of particles created outside region \"Target\" :" << G4endl;
  PrintNParticles(particlesCreatedInWorld);

  G4cout << "_______________________" << G4endl;
  std::map<const G4ParticleDefinition*, int>&  particlesCreatedInTarget = TrackingAct->GetNParticlesCreatedInTarget();
  G4cout << "Number and type of particles created in region \"Target\" :" << G4endl;
  PrintNParticles(particlesCreatedInTarget);
}

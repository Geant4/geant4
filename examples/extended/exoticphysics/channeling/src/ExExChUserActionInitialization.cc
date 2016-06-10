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

#ifdef G4MULTITHREADED
#include "ExExChUserActionInitialization.hh"
#include "ExExChPrimaryGeneratorAction.hh"
#include "ExExChTrackingAction.hh"
#include "ExExChStackingAction.hh"
#include "ExExChEventAction.hh"
#include "ExExChRunAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExExChUserActionInitialization::ExExChUserActionInitialization() {
    //GPS is special: it is thread-private,
    //but underlying parameters and objects governing
    //the distributions are shared among threads.
    //We need to be sure that the defaults are set only once by master
    masterGPS = new G4GeneralParticleSource();
    masterGPS->SetParticleDefinition(G4ParticleTable::
                        GetParticleTable()->FindParticle("proton"));

    // Position distribution
    G4SPSPosDistribution *vPosDist =
        masterGPS->GetCurrentSource()->GetPosDist();
    vPosDist->SetPosDisType("Beam");
    vPosDist->SetPosDisShape("Circle");
    vPosDist->SetCentreCoords(G4ThreeVector(0.,0.,-1.05 * CLHEP::meter));
    vPosDist->SetBeamSigmaInR(0.0 * CLHEP::mm);

    // Angular distribution
    G4SPSAngDistribution *vAngDist =
        masterGPS->GetCurrentSource()->GetAngDist();
    vAngDist->DefineAngRefAxes("angref1",G4ThreeVector(1.,0.,0));
    vAngDist->DefineAngRefAxes("angref2",G4ThreeVector(0.,-1.,0));
    vAngDist->SetAngDistType("beam2d");
    vAngDist->SetBeamSigmaInAngX(13.36E-6 * CLHEP::rad);
    vAngDist->SetBeamSigmaInAngY(11.25E-6 * CLHEP::rad);

    // Energy distribution
    G4SPSEneDistribution *vEneDist =
        masterGPS->GetCurrentSource()->GetEneDist();
    vEneDist->SetEnergyDisType("Mono");
    vEneDist->SetMonoEnergy(400. * CLHEP::GeV);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExExChUserActionInitialization::~ExExChUserActionInitialization() {
  delete masterGPS;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChUserActionInitialization::Build() const {
    SetUserAction(new ExExChPrimaryGeneratorAction());
    SetUserAction(new ExExChEventAction());
    SetUserAction(new ExExChStackingAction());
    SetUserAction(new ExExChTrackingAction());
    SetUserAction(new ExExChRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChUserActionInitialization::BuildForMaster() const {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#endif

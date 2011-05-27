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
#include "G4RunManager.hh"
#include "DummyDetectorConstruction.hh"
#include "DummyPhysicsList.hh"
#include "DummyPrimaryGeneratorAction.hh"

#include "TFile.h"
#include "TH1F.h"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"

#include "G4HadronCrossSections.hh"
#include "G4HadronFissionProcess.hh"
#include "G4LFission.hh"


int main(int argc,char** argv)
{
  if (argc != 3) {
    G4cout << " Projectile momentum (GeV/c) and event number arguments expected " << G4endl;
    return 0;
  }
  G4double projectileMomentum = atof(argv[1]);
  G4int Nevents = atoi(argv[2]);

  G4cout << " neutron momentum = " << projectileMomentum << " GeV/c " << G4endl;
  G4cout << " N = " << Nevents << G4endl;

  // RunManager construction (only necessary for particle table lookup)
  G4RunManager* runManager = new G4RunManager;

  // mandatory user initialization classes
  runManager->SetUserInitialization(new DummyDetectorConstruction);
  runManager->SetUserInitialization(new DummyPhysicsList);

  // initialize Geant4 kernel
  runManager->Initialize();

  // mandatory user action class
  runManager->SetUserAction(new DummyPrimaryGeneratorAction);

  /////////////////////////////////////////////////////////////
  //                                                         //
  //     Set up histograms                                   //
  //                                                         //
  /////////////////////////////////////////////////////////////

  TH1F* cost =  new TH1F("cost", "cos(theta) of proton",110,0.0, 1.1);
  TH1F* pang =  new TH1F("pang", "theta",180,0.0, 180.0);
  TH1F* pkin =  new TH1F("pkin", "proton KE",100,0.0,10.0);

  /////////////////////////////////////////////////////////////
  //                                                         //
  //     Define material                                     // 
  //                                                         //
  /////////////////////////////////////////////////////////////

  std::vector<G4Material*> materials;
  materials.push_back(new G4Material("Th", 90.0, 232.0*g/mole, 11.7*g/cm3) );
  materials.push_back(new G4Material("U235", 92.0, 235.0*g/mole, 18.95*g/cm3) );

  /////////////////////////////////////////////////////////////
  //                                                         //
  //     Set up capture process                              //
  //                                                         //
  /////////////////////////////////////////////////////////////

  G4ParticleDefinition* part = G4Neutron::Neutron();
  G4HadronFissionProcess* fissionProcess = new G4HadronFissionProcess();
  fissionProcess->RegisterMe(new G4LFission() );

  G4ThreeVector aMomentum = G4ThreeVector(0.,0.,projectileMomentum*GeV);
  G4ThreeVector aPosition(0., 0., 0.);
  G4double aTime = 0.;
  G4ForceCondition cond = Forced;
  G4StepPoint* aPoint = new G4StepPoint();
  aPoint->SetPosition(aPosition);
  G4Step step;
  G4VParticleChange* aChange;

  // Define the track

  G4DynamicParticle dParticle(part,aMomentum);
  G4cout << part->GetParticleName() << " KE = " 
         << dParticle.GetKineticEnergy()/GeV << G4endl;
  G4Track* track = new G4Track(&dParticle,aTime,aPosition);
  step.SetTrack(track);
  TFile* hfile = 0;
  G4String matName = "H";

  for (G4int i = 0; i < materials.size(); i++) {
    matName = materials[i]->GetName();
    if (matName == "Th") {
      hfile = new TFile("Th.root","new");
    } else if (matName == "U235") {
      hfile = new TFile("U235.root","new");
    } else {
      return 0;
    }

    aPoint->SetMaterial(materials[i]);
    step.SetPreStepPoint(aPoint); 
    track->SetStep(&step);

    /////////////////////////////////////////////////////////////
    //                                                         //
    //   Set up event loop                                     //
    //                                                         //
    /////////////////////////////////////////////////////////////

    G4Track* sec = 0;
    G4int nSec = 0;
    G4int nSecTot = 0;

    G4double npx = 0.;
    G4double npy = 0.;
    G4double npz = 0.;
    G4double nke = 0.;
    G4int ntot = 0;

    G4double gpx = 0.;
    G4double gpy = 0.;
    G4double gpz = 0.;
    G4double gke = 0.;
    G4int gtot = 0;

    for (G4int loop = 0; loop < Nevents; loop++) {
      aChange = fissionProcess->PostStepDoIt(*track,step);
      nSec = aChange->GetNumberOfSecondaries();
      nSecTot += nSec;
      for (G4int i = 0; i < nSec; i++) {
        sec = aChange->GetSecondary(i);
        if (sec->GetDefinition()->GetParticleName() == "neutron") {
          npx += sec->GetMomentumDirection().x();
          npy += sec->GetMomentumDirection().y();
          npz += sec->GetMomentumDirection().z();
          nke += sec->GetKineticEnergy();
          ntot++;
	} else if (sec->GetDefinition()->GetParticleName() == "gamma") {
          gpx += sec->GetMomentumDirection().x();
          gpy += sec->GetMomentumDirection().y();
          gpz += sec->GetMomentumDirection().z();
          gke += sec->GetKineticEnergy();
          gtot++;
        }
      }
    }
 
    G4cout << G4endl << " For a " << matName << " target, the total number of gammas = " << gtot << G4endl;
    G4cout << G4double(gtot)/G4double(Nevents) << " gammas per event " << G4endl;
    G4cout << " mean kinetic energy = " << gke/G4double(Nevents)
           << " , mean px/p = " << gpx/G4double(Nevents) 
           << " , mean py/p = " << gpy/G4double(Nevents)
           << " , mean pz/p = " << gpz/G4double(Nevents) << G4endl;
 
    G4cout << " The total number of neutrons = " << ntot << G4endl;
    G4cout << G4double(ntot)/G4double(Nevents) << " neutrons per event " << G4endl;
    G4cout << " mean kinetic energy = " << nke/G4double(Nevents)
           << " , mean px/p = " << npx/G4double(Nevents) 
           << " , mean py/p = " << npy/G4double(Nevents)
           << " , mean pz/p = " << npz/G4double(Nevents) << G4endl;  
  }

  /////////////////////////////////////////////////////////////
  //                                                         //
  //     Save histograms                                     //
  //                                                         //
  /////////////////////////////////////////////////////////////

  hfile->Write();
  hfile->Close();

  delete runManager;
  return 0;
}

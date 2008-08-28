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
#include "G4HadronElasticProcess.hh"
#include "G4LElastic.hh"


int main(int argc,char** argv)
{
  if (argc != 3) {
    G4cout << " Projectile momentum (GeV/c) and event number arguments expected " << G4endl;
    return 0;
  }
  G4double projectileMomentum = atof(argv[1]);
  G4int Nevents = atoi(argv[2]);

  G4cout << " proton momentum = " << projectileMomentum << " GeV/c " << G4endl;
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
  materials.push_back(new G4Material("H", 1.0, 1.008*g/mole, 0.071*g/cm3) );
  materials.push_back(new G4Material("Be", 4.0, 9.012*g/mole, 1.848*g/cm3) );
  materials.push_back(new G4Material("C", 6.0, 12.01*g/mole, 2.265*g/cm3) );
  materials.push_back(new G4Material("Fe", 26.0, 55.845*g/mole, 7.87*g/cm3) );
  materials.push_back(new G4Material("U235", 92.0, 235.0*g/mole, 18.95*g/cm3) );

  /////////////////////////////////////////////////////////////
  //                                                         //
  //     Set up elastic process                              //
  //                                                         //
  /////////////////////////////////////////////////////////////

  G4ParticleDefinition* part = G4Proton::Proton();
  G4HadronElasticProcess* elasticProcess = new G4HadronElasticProcess();
  elasticProcess->RegisterMe(new G4LElastic() );

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
    if (matName == "H") {
      hfile = new TFile("H.root","new");
    } else if (matName == "Be") {
      hfile = new TFile("Be.root","new");
    } else if (matName == "C") {
      hfile = new TFile("C.root","new");
    } else if (matName == "Fe") {
      hfile = new TFile("Fe.root","new");
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
    G4String secName;
    G4ThreeVector secMom;

    G4double pxg = 0;
    G4double pyg = 0;
    G4double pzg = 0;

    G4double pnew = 0;
    G4double pxnew = 0;
    G4double pynew = 0;
    G4double pznew = 0;

    G4double ekg = 0;
    G4double eknew4 = 0;
    G4double pxnew4 = 0;
    G4double pynew4 = 0; 
    G4double pznew4 = 0;

    G4int nSecTot = 0;
    G4double costTot = 0.;
    G4double KETot = 0.;

    for (G4int loop = 0; loop < Nevents; loop++) {
      aChange = elasticProcess->PostStepDoIt(*track,step);

      // Get end of step information which contains momentum of
      // scattered proto
      G4Step* updatedStep = aChange->UpdateStepForPostStep(&step);
      G4StepPoint* endOfStep = updatedStep->GetPostStepPoint();
      costTot += endOfStep->GetMomentumDirection().z();
      KETot += endOfStep->GetKineticEnergy();
      nSecTot += aChange->GetNumberOfSecondaries();
    }
 
    G4cout << G4endl << " For a " << matName << " target, the mean KE of the scattered proton = " 
           << KETot/G4double(Nevents) << G4endl;
    G4cout << " mean value of cos(theta) = " << costTot/G4double(Nevents) << G4endl;
    G4cout << " mean number of secondaries = " << G4double(nSecTot)/G4double(Nevents) << G4endl;
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

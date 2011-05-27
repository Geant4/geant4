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
//   File:    HEInelastic.cc - updated process-level test of HEP models
//   Author:  D.H. Wright (SLAC)
//   Date:    29 Aug 2008
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

#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4LambdaInelasticProcess.hh"
#include "G4SigmaPlusInelasticProcess.hh"
#include "G4SigmaMinusInelasticProcess.hh"
#include "G4XiZeroInelasticProcess.hh"
#include "G4XiMinusInelasticProcess.hh"
#include "G4OmegaMinusInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4AntiLambdaInelasticProcess.hh"
#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4AntiXiZeroInelasticProcess.hh"
#include "G4AntiXiMinusInelasticProcess.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"

#include "G4HEProtonInelastic.hh"
#include "G4HENeutronInelastic.hh"
#include "G4HELambdaInelastic.hh"
#include "G4HESigmaPlusInelastic.hh"
#include "G4HESigmaMinusInelastic.hh"
#include "G4HEXiZeroInelastic.hh"
#include "G4HEXiMinusInelastic.hh"
#include "G4HEOmegaMinusInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4HEAntiNeutronInelastic.hh"
#include "G4HEAntiLambdaInelastic.hh"
#include "G4HEAntiSigmaPlusInelastic.hh"
#include "G4HEAntiSigmaMinusInelastic.hh"
#include "G4HEAntiXiZeroInelastic.hh"
#include "G4HEAntiXiMinusInelastic.hh"
#include "G4HEAntiOmegaMinusInelastic.hh"
#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"


int main(int argc,char** argv)
{
  if (argc != 3) {
    G4cout << " Projectile momentum (GeV/c) and event number arguments expected " << G4endl;
    return 0;
  }
  G4double projectileMomentum = atof(argv[1]);
  G4int Nevents = atoi(argv[2]);
  G4double nevt = Nevents;

  G4cout << " initial momentum = " << projectileMomentum << " GeV/c " << G4endl;
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
  //  materials.push_back(new G4Material("H", 1.0, 1.008*g/mole, 0.071*g/cm3) );
  materials.push_back(new G4Material("C", 6.0, 12.01*g/mole, 2.265*g/cm3) );
  materials.push_back(new G4Material("Fe", 26.0, 55.845*g/mole, 7.87*g/cm3) );
  materials.push_back(new G4Material("Pb", 82.0, 207.2*g/mole, 11.34*g/cm3) );

  /////////////////////////////////////////////////////////////
  //                                                         //
  //     Set up elastic process                              //
  //                                                         //
  /////////////////////////////////////////////////////////////

  G4HadronInelasticProcess* inelasticProcess = 0;
  std::vector<G4ParticleDefinition*> particles;
  std::vector<G4HadronInelasticProcess*> processes;

  particles.push_back(G4Proton::Proton() );
  inelasticProcess = new G4ProtonInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEProtonInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4Neutron::Neutron() );
  inelasticProcess = new G4NeutronInelasticProcess();
  inelasticProcess->RegisterMe(new G4HENeutronInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4Lambda::Lambda() );
  inelasticProcess = new G4LambdaInelasticProcess();
  inelasticProcess->RegisterMe(new G4HELambdaInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4SigmaPlus::SigmaPlus() );
  inelasticProcess = new G4SigmaPlusInelasticProcess();
  inelasticProcess->RegisterMe(new G4HESigmaPlusInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4SigmaMinus::SigmaMinus() );
  inelasticProcess = new G4SigmaMinusInelasticProcess();
  inelasticProcess->RegisterMe(new G4HESigmaMinusInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4XiZero::XiZero() );
  inelasticProcess = new G4XiZeroInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEXiZeroInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4XiMinus::XiMinus() );
  inelasticProcess = new G4XiMinusInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEXiMinusInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4OmegaMinus::OmegaMinus() );
  inelasticProcess = new G4OmegaMinusInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEOmegaMinusInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4AntiProton::AntiProton() );
  inelasticProcess = new G4AntiProtonInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEAntiProtonInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4AntiNeutron::AntiNeutron() );
  inelasticProcess = new G4AntiNeutronInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEAntiNeutronInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4AntiLambda::AntiLambda() );
  inelasticProcess = new G4AntiLambdaInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEAntiLambdaInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4AntiSigmaPlus::AntiSigmaPlus() );
  inelasticProcess = new G4AntiSigmaPlusInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEAntiSigmaPlusInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4AntiSigmaMinus::AntiSigmaMinus() );
  inelasticProcess = new G4AntiSigmaMinusInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEAntiSigmaMinusInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4AntiXiZero::AntiXiZero() );
  inelasticProcess = new G4AntiXiZeroInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEAntiXiZeroInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4AntiXiMinus::AntiXiMinus() );
  inelasticProcess = new G4AntiXiMinusInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEAntiXiMinusInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4AntiOmegaMinus::AntiOmegaMinus() );
  inelasticProcess = new G4AntiOmegaMinusInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEAntiOmegaMinusInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4PionPlus::PionPlus() );
  inelasticProcess = new G4PionPlusInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEPionPlusInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4PionMinus::PionMinus() );
  inelasticProcess = new G4PionMinusInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEPionMinusInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4KaonPlus::KaonPlus() );
  inelasticProcess = new G4KaonPlusInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEKaonPlusInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4KaonMinus::KaonMinus() );
  inelasticProcess = new G4KaonMinusInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEKaonMinusInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4KaonZeroLong::KaonZeroLong() );
  inelasticProcess = new G4KaonZeroLInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEKaonZeroInelastic() );
  processes.push_back(inelasticProcess);

  particles.push_back(G4KaonZeroShort::KaonZeroShort() );
  inelasticProcess = new G4KaonZeroSInelasticProcess();
  inelasticProcess->RegisterMe(new G4HEKaonZeroInelastic() );
  processes.push_back(inelasticProcess);

  G4ThreeVector aMomentum = G4ThreeVector(0.,0.,projectileMomentum*GeV);
  G4ThreeVector aPosition(0., 0., 0.);
  G4double aTime = 0.;
  G4StepPoint* aPoint = new G4StepPoint();
  aPoint->SetPosition(aPosition);
  G4Step step;
  G4VParticleChange* aChange;

  /*
  TFile* hfile = 0;
  G4String matName = "H";
  for (G4int i = 0; i < materials.size(); i++) {
    matName = materials[i]->GetName();
    if (matName == "H") {
      hfile = new TFile("H.root","new");
    } else if (matName == "C") {
      hfile = new TFile("C.root","new");
    } else if (matName == "Fe") {
      hfile = new TFile("Fe.root","new");
    } else if (matName == "Pb") {
      hfile = new TFile("Pb.root","new");
    } else {
      return 0;
    }
  }
  */

  G4Track* track = 0;
  G4Track* sec = 0;
  G4String secName;
  G4ParticleDefinition* pDef = 0;

  // Loop over particle types
  for (G4int i = 0; i < particles.size(); i++) {
    G4DynamicParticle dParticle(particles[i],aMomentum);
    G4cout << G4endl << projectileMomentum << " GeV " << particles[i]->GetParticleName() << G4endl;

    // Loop over materials
    for (G4int j = 0; j < materials.size(); j++) {
      G4cout << " in " << materials[j]->GetName() << G4endl;
      // Define track and step
      track = new G4Track(&dParticle,aTime,aPosition);
      step.SetTrack(track);
      aPoint->SetMaterial(materials[j]);
      step.SetPreStepPoint(aPoint); 
      track->SetStep(&step);

      // Set up event loop

      G4double KE = 0;

      G4int nprot = 0;
      G4int nneut = 0;
      G4int npip = 0;
      G4int npim = 0;
      G4int npiz = 0;
      G4int nkp = 0;
      G4int nkm = 0;
      G4int nkz = 0;
      G4int nlam = 0;
      G4int npbar = 0;

      G4double sumKEp = 0;
      G4double sumKEn = 0.;
      G4double sumKEpip = 0.;
      G4double sumKEpim = 0.;
      G4double sumKEpiz = 0.;
      G4double sumKEkp = 0.;
      G4double sumKEkm = 0.;
      G4double sumKEkz = 0.;
      G4double sumKEl = 0.;
      G4double sumKEpbar = 0.;

      G4int nSec = 0;
      G4int nSecTot = 0;
      G4double costTot = 0.;
      G4double KETot = 0.;

      for (G4int loop = 0; loop < Nevents; loop++) {
        aChange = inelasticProcess->PostStepDoIt(*track,step);
        nSec = aChange->GetNumberOfSecondaries();
        nSecTot += nSec;
        for (G4int k = 0; k < nSec; k++) {
          sec = aChange->GetSecondary(k);
          KE = sec->GetKineticEnergy()/GeV;
          pDef = sec->GetDefinition();
          if (pDef == G4Proton::Proton()) {
            nprot++;
            sumKEp += KE;
          } else if (pDef == G4Neutron::Neutron()) {
            nneut++;
            sumKEn += KE;
          } else if (pDef == G4Lambda::Lambda()) {
            nlam++;
            sumKEl += KE;
          } else if (pDef == G4AntiProton::AntiProton()) {
            npbar++;
            sumKEpbar += KE;
          } else if (pDef == G4PionPlus::PionPlus()) {
            npip++;
            sumKEpip += KE;
          } else if (pDef == G4PionMinus::PionMinus()) {
            npim++;
            sumKEpim += KE;
          } else if (pDef == G4PionZero::PionZero()) {
            npiz++;
            sumKEpiz += KE;
          } else if (pDef == G4KaonPlus::KaonPlus()) {
            nkp++;
            sumKEkp += KE;
          } else if (pDef == G4KaonMinus::KaonMinus()) {
            nkm++;
            sumKEkm += KE;
          } else if (pDef == G4KaonZero::KaonZero() || 
                     pDef == G4KaonZeroLong::KaonZeroLong() ||
                     pDef == G4KaonZeroShort::KaonZeroShort() ) {
            nkz++;
            sumKEkz += KE;
          }
        }
      }  // event loop
      G4cout << G4double(nprot)/nevt << " protons per event with mean KE " 
             << sumKEp/nevt << " GeV " << G4endl;
      G4cout << G4double(nneut)/nevt << " neutrons per event with mean KE " 
             << sumKEn/nevt << " GeV " << G4endl;
      G4cout << G4double(nlam)/nevt << " lambdas per event with mean KE " 
             << sumKEl/nevt << " GeV " << G4endl;
      G4cout << G4double(npbar)/nevt << " anti-protons per event with mean KE " 
             << sumKEpbar/nevt << " GeV " << G4endl;
      G4cout << G4double(npip)/nevt << " pi plus per event with mean KE " 
             << sumKEpip/nevt << " GeV " << G4endl;
      G4cout << G4double(npim)/nevt << " pi minus per event with mean KE " 
             << sumKEpim/nevt << " GeV " << G4endl;
      G4cout << G4double(npiz)/nevt << " pi zero per event with mean KE " 
             << sumKEpiz/nevt << " GeV " << G4endl;
      G4cout << G4double(nkp)/nevt << " K plus per event with mean KE " 
             << sumKEkp/nevt << " GeV " << G4endl;
      G4cout << G4double(nkm)/nevt << " K minus per event with mean KE " 
             << sumKEkm/nevt << " GeV " << G4endl;
      G4cout << G4double(nkz)/nevt << " K0 per event with mean KE " 
             << sumKEkz/nevt << " GeV " << G4endl;
    }  // material loop
  } // particle loop
 
  /////////////////////////////////////////////////////////////
  //                                                         //
  //     Save histograms                                     //
  //                                                         //
  /////////////////////////////////////////////////////////////

  //  hfile->Write();
  //  hfile->Close();

  delete runManager;
  return 0;
}

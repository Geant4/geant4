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
//---------------------------------------------------------------------------
//
// ClassName:   G4NeutrinoPhysics
//
// Author: 2023 V. Ivanchenko extracted from G4EmExtraPhysics
//
// Modified:
//
//
///////////////////////////////////////////////////////////////

#include "G4NeutrinoPhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4Electron.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4NeutrinoTau.hh"

#include "G4NeutrinoElectronProcess.hh"
#include "G4NeutrinoElectronTotXsc.hh"
#include "G4NeutrinoElectronCcModel.hh"
#include "G4NeutrinoElectronNcModel.hh"

#include "G4MuNeutrinoNucleusProcess.hh"
#include "G4TauNeutrinoNucleusProcess.hh"
#include "G4ElNeutrinoNucleusProcess.hh"
#include "G4NuVacOscProcess.hh"

#include "G4MuNeutrinoNucleusTotXsc.hh"
#include "G4TauNeutrinoNucleusTotXsc.hh"
#include "G4ElNeutrinoNucleusTotXsc.hh"

#include "G4NuMuNucleusCcModel.hh"
#include "G4NuMuNucleusNcModel.hh"
#include "G4ANuMuNucleusCcModel.hh"
#include "G4ANuMuNucleusNcModel.hh"

#include "G4NuTauNucleusCcModel.hh"
#include "G4NuTauNucleusNcModel.hh"
#include "G4ANuTauNucleusCcModel.hh"
#include "G4ANuTauNucleusNcModel.hh"

#include "G4NuElNucleusCcModel.hh"
#include "G4NuElNucleusNcModel.hh"
#include "G4ANuElNucleusCcModel.hh"
#include "G4ANuElNucleusNcModel.hh"
 
// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4NeutrinoPhysics);

//////////////////////////////////////

G4NeutrinoPhysics::G4NeutrinoPhysics(G4int ver): 
  G4VPhysicsConstructor("NeutrinoPhys"),
  verbose(ver)
{
  theMessenger = new G4NeutrinoPhysicsMessenger(this);
  if(verbose > 1) G4cout << "### G4NeutrinoPhysics" << G4endl;
}

G4NeutrinoPhysics::~G4NeutrinoPhysics()
{
  delete theMessenger;
}

void G4NeutrinoPhysics::NuETotXscActivated(G4bool val)
{
  fNuETotXscActivated = val;
}

void G4NeutrinoPhysics::SetNuOscillation(G4bool val)
{
  fNuOscillation = val;
}

void G4NeutrinoPhysics::SetNuEleCcBias(G4double bf)
{
  if(bf > 0.0) fNuEleCcBias = bf;
}

void G4NeutrinoPhysics::SetNuEleNcBias(G4double bf)
{
  if(bf > 0.0) fNuEleNcBias = bf;
}

void G4NeutrinoPhysics::SetNuNucleusBias(G4double bf)
{
  if(bf > 0.0) fNuNucleusBias = bf;
}

void G4NeutrinoPhysics::SetNuOscDistanceBias(G4double bf)
{
  if(bf > 0.0) fNuOscDistanceBias = bf;
}

void G4NeutrinoPhysics::SetNuDetectorName(const G4String& dn)
{
  fNuDetectorName = dn;
}

void G4NeutrinoPhysics::SetNuOscDistanceName(const G4String& dn)
{
  fNuOscDistanceName = dn;
}

/////////////////////////////////////////////////

void G4NeutrinoPhysics::ConstructParticle()
{
  G4Electron::Electron();
  G4AntiNeutrinoE::AntiNeutrinoE();
  G4NeutrinoE::NeutrinoE();
  G4AntiNeutrinoMu::AntiNeutrinoMu();
  G4NeutrinoMu::NeutrinoMu();
  G4AntiNeutrinoTau::AntiNeutrinoTau();
  G4NeutrinoTau::NeutrinoTau();
}

void G4NeutrinoPhysics::ConstructProcess()
{
  const G4ParticleDefinition* p[6] = {
    G4AntiNeutrinoE::AntiNeutrinoE(),
    G4NeutrinoE::NeutrinoE(),
    G4AntiNeutrinoMu::AntiNeutrinoMu(),
    G4NeutrinoMu::NeutrinoMu(),
    G4AntiNeutrinoTau::AntiNeutrinoTau(),
    G4NeutrinoTau::NeutrinoTau()
  };

  // neutrino vacuum oscillation process
  if (fNuOscillation) {
    auto theNuVacOscProcess = new G4NuVacOscProcess(fNuOscDistanceName);
    theNuVacOscProcess->SetBiasingFactor(fNuOscDistanceBias);
    

    for (G4int i=0; i<6; ++i) {
      p[i]->GetProcessManager()->AddDiscreteProcess(theNuVacOscProcess);
    }
  }

  // neutrino-electron process
  auto theNuEleProcess = new G4NeutrinoElectronProcess(fNuDetectorName);
  G4NeutrinoElectronTotXsc* theNuEleTotXsc = new G4NeutrinoElectronTotXsc();

  if (fNuETotXscActivated) {
    G4double bftot = std::max(fNuEleCcBias, fNuEleNcBias);
    theNuEleProcess->SetBiasingFactor(bftot);
  }
  else {
    theNuEleProcess->SetBiasingFactors(fNuEleCcBias, fNuEleNcBias);
    theNuEleTotXsc->SetBiasingFactors(fNuEleCcBias, fNuEleNcBias);
  }
  theNuEleProcess->AddDataSet(theNuEleTotXsc);

  G4NeutrinoElectronCcModel* ccModel = new G4NeutrinoElectronCcModel();
  G4NeutrinoElectronNcModel* ncModel = new G4NeutrinoElectronNcModel();
  theNuEleProcess->RegisterMe(ccModel);
  theNuEleProcess->RegisterMe(ncModel);

  for (G4int i=0; i<6; ++i) {
    p[i]->GetProcessManager()->AddDiscreteProcess(theNuEleProcess);
  }

  // nu_mu nucleus interactions
  auto theNuMuNucleusProcess = new G4MuNeutrinoNucleusProcess(fNuDetectorName);
  auto theNuMuNucleusTotXsc = new G4MuNeutrinoNucleusTotXsc();
    
  if (fNuETotXscActivated) {
    theNuMuNucleusProcess->SetBiasingFactor(fNuNucleusBias);
  }
  theNuMuNucleusProcess->AddDataSet(theNuMuNucleusTotXsc);

  G4NuMuNucleusCcModel* numunuclcc = new G4NuMuNucleusCcModel();
  G4NuMuNucleusNcModel* numunuclnc = new G4NuMuNucleusNcModel();
  G4ANuMuNucleusCcModel* anumunuclcc = new G4ANuMuNucleusCcModel();
  G4ANuMuNucleusNcModel* anumunuclnc = new G4ANuMuNucleusNcModel();
    
  theNuMuNucleusProcess->RegisterMe(numunuclcc);
  theNuMuNucleusProcess->RegisterMe(numunuclnc);
  theNuMuNucleusProcess->RegisterMe(anumunuclcc);
  theNuMuNucleusProcess->RegisterMe(anumunuclnc);

  for (G4int i=2; i<=3; ++i) {
    p[i]->GetProcessManager()->AddDiscreteProcess(theNuMuNucleusProcess);
  }

  // nu_tau nucleus interactions
  auto theNuTauNucleusProcess = new G4TauNeutrinoNucleusProcess(fNuDetectorName);
  auto theNuTauNucleusTotXsc = new G4TauNeutrinoNucleusTotXsc();
    
  if(fNuETotXscActivated) {
    theNuTauNucleusProcess->SetBiasingFactor(fNuNucleusBias);
  }
  theNuTauNucleusProcess->AddDataSet(theNuTauNucleusTotXsc);

  G4NuTauNucleusCcModel* nutaunuclcc = new G4NuTauNucleusCcModel();
  G4NuTauNucleusNcModel* nutaunuclnc = new G4NuTauNucleusNcModel();
  G4ANuTauNucleusCcModel* anutaunuclcc = new G4ANuTauNucleusCcModel();
  G4ANuTauNucleusNcModel* anutaunuclnc = new G4ANuTauNucleusNcModel();
    
  theNuTauNucleusProcess->RegisterMe(nutaunuclcc);
  theNuTauNucleusProcess->RegisterMe(nutaunuclnc);
  theNuTauNucleusProcess->RegisterMe(anutaunuclcc);
  theNuTauNucleusProcess->RegisterMe(anutaunuclnc);

  for (G4int i=4; i<=5; ++i) {
    p[i]->GetProcessManager()->AddDiscreteProcess(theNuMuNucleusProcess);
  }

  // nu_e nucleus interactions
  auto theNuElNucleusProcess = new G4ElNeutrinoNucleusProcess(fNuDetectorName);
  auto theNuElNucleusTotXsc = new G4ElNeutrinoNucleusTotXsc();
    
  if (fNuETotXscActivated) {
    theNuElNucleusProcess->SetBiasingFactor(fNuNucleusBias);
  }
  theNuElNucleusProcess->AddDataSet(theNuElNucleusTotXsc);

  G4NuElNucleusCcModel* nuelnuclcc = new G4NuElNucleusCcModel();
  G4NuElNucleusNcModel* nuelnuclnc = new G4NuElNucleusNcModel();
  G4ANuElNucleusCcModel* anuelnuclcc = new G4ANuElNucleusCcModel();
  G4ANuElNucleusNcModel* anuelnuclnc = new G4ANuElNucleusNcModel();
    
  theNuElNucleusProcess->RegisterMe(nuelnuclcc);
  theNuElNucleusProcess->RegisterMe(nuelnuclnc);
  theNuElNucleusProcess->RegisterMe(anuelnuclcc);
  theNuElNucleusProcess->RegisterMe(anuelnuclnc);

  for (G4int i=0; i<=1; ++i) {
    p[i]->GetProcessManager()->AddDiscreteProcess(theNuElNucleusProcess);
  }
}

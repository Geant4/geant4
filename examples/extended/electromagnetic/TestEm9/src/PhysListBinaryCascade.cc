//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: PhysListBinaryCascade.cc,v 1.4 2003/12/05 11:22:45 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $

#include "PhysListBinaryCascade.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

#include "G4BinaryCascade.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListBinaryCascade::PhysListBinaryCascade(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListBinaryCascade::~PhysListBinaryCascade()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListBinaryCascade::ConstructProcess()
{

  // Binary Cascade
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* pmanager = 0;

  G4BinaryCascade* theBC = 0;

  // proton
  particle = G4Proton::Proton();
  pmanager = particle->GetProcessManager();
  theBC = new G4BinaryCascade();
  theIPproton.RegisterMe(theBC);
  theIPproton.AddDataSet(&thePXSec);
  pmanager->AddDiscreteProcess(&theIPproton);

  // neutron
  particle = G4Neutron::Neutron();
  pmanager = particle->GetProcessManager();
  theBC = new G4BinaryCascade();
  theIPneutron.RegisterMe(theBC);
  theIPneutron.AddDataSet(&theNXSec);
  pmanager->AddDiscreteProcess(&theIPneutron);
  // fission
  G4HadronFissionProcess* theFissionProcess = new G4HadronFissionProcess;
  G4LFission* theFissionModel = new G4LFission;
  theFissionProcess->RegisterMe(theFissionModel);
  pmanager->AddDiscreteProcess(theFissionProcess);
  // capture
  G4HadronCaptureProcess* theCaptureProcess = new G4HadronCaptureProcess;
  G4LCapture* theCaptureModel = new G4LCapture;
  theCaptureProcess->RegisterMe(theCaptureModel);
  pmanager->AddDiscreteProcess(theCaptureProcess);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


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
// $Id: PhysListBinaryCascade.cc,v 1.7 2006/06/29 16:58:19 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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


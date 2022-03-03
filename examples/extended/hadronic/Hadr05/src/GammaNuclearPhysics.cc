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
/// \file GammaNuclearPhysics.cc
/// \brief Implementation of the GammaNuclearPhysics class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "GammaNuclearPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

// Processes

#include "G4HadronInelasticProcess.hh"
#include "G4LowEGammaNuclearModel.hh"
#include "G4CascadeInterface.hh"

#include "G4PhotoNuclearCrossSection.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GammaNuclearPhysics::GammaNuclearPhysics(const G4String& name)
:  G4VPhysicsConstructor(name)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GammaNuclearPhysics::~GammaNuclearPhysics()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GammaNuclearPhysics::ConstructProcess()
{
   G4HadronInelasticProcess* process
       = new G4HadronInelasticProcess("photonNuclear", G4Gamma::Definition());
   process->AddDataSet( new G4PhotoNuclearCrossSection );   

   // to not register a model, set Emax=0; eg. Emax1 = 0.
   const G4double Emax1 = 200*MeV, Emax2 = 10*GeV;

   if (Emax1 > 0.) {  // model 1
     G4LowEGammaNuclearModel* model1 = new G4LowEGammaNuclearModel();
     model1->SetMaxEnergy(Emax1);
     process->RegisterMe(model1);
   }
   
   if (Emax2 > 0.) {  // model 2   
     G4CascadeInterface* model2 = new G4CascadeInterface();
     G4double Emin2 = std::max(Emax1-1*MeV, 0.);
     model2->SetMinEnergy(Emin2);   
     model2->SetMaxEnergy(Emax2);
     process->RegisterMe(model2);
   }
   
   G4ProcessManager* pManager = G4Gamma::Gamma()->GetProcessManager();
   pManager->AddDiscreteProcess(process);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

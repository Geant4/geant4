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
/// \file hadronic/Hadr02/src/UrQMDNeutronBuilder.cc
/// \brief Implementation of the UrQMDNeutronBuilder class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   UrQMDNeutronBuilder
//
// Author: 2012 Andrea Dotti
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifdef G4_USE_URQMD
#include "UrQMDNeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4HadronicParameters.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UrQMDNeutronBuilder::UrQMDNeutronBuilder() 
{
  fMin = 0*MeV;
  fMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  fModel = new G4UrQMD1_3Model();
  captureModel = new G4NeutronRadCapture();
  fissionModel = new G4LFission();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UrQMDNeutronBuilder::Build(G4NeutronInelasticProcess * aP)
{
  fModel->SetMinEnergy(fMin);
  fModel->SetMaxEnergy(fMax);
  aP->RegisterMe(fModel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UrQMDNeutronBuilder::~UrQMDNeutronBuilder() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UrQMDNeutronBuilder::Build(G4HadronElasticProcess * )
{}

void UrQMDNeutronBuilder::Build(G4HadronFissionProcess* aP)
{
  fissionModel->SetMinEnergy(0.0);
  fissionModel->SetMaxEnergy(20.0*TeV);
  aP->RegisterMe(fissionModel);
}

void UrQMDNeutronBuilder::Build(G4HadronCaptureProcess* aP)
{
  aP->RegisterMe(captureModel);
}

#endif //G4_USE_URQMD

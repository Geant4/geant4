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
/// \file hadronic/Hadr02/src/UrQMDAntiBarionBuilder.cc
/// \brief Implementation of the UrQMDAntiBarionBuilder class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   UrQMDAntiBarionBuilder
//
// Author: 2012 Andrea Dotti
//
// Modified:
//----------------------------------------------------------------------------
//
#ifdef G4_USE_URQMD
#include "UrQMDAntiBarionBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"  // For anti-ions
#include "G4CrossSectionInelastic.hh"
#include "G4HadronicParameters.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UrQMDAntiBarionBuilder::UrQMDAntiBarionBuilder() 
{
  //Set-up UrQMD model
  fMin =   0.0*MeV;
  fMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  fModel = new G4UrQMD1_3Model();
  fModel->SetMinEnergy( fMin );
  fModel->SetMaxEnergy( fMax );
  fAntiNucleonXS=new G4ComponentAntiNuclNuclearXS();

  fAntiNucleonData = 
    new G4CrossSectionInelastic(fAntiNucleonXS);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UrQMDAntiBarionBuilder::~UrQMDAntiBarionBuilder() 
{
  delete fAntiNucleonXS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UrQMDAntiBarionBuilder::Build(G4HadronElasticProcess * ) 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UrQMDAntiBarionBuilder::Build(G4AntiProtonInelasticProcess * aP)
{
  fModel->SetMinEnergy(fMin);
  fModel->SetMaxEnergy(fMax);
  aP->AddDataSet(fAntiNucleonData);
  aP->RegisterMe(fModel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UrQMDAntiBarionBuilder::Build(G4AntiNeutronInelasticProcess * aP)
{
  fModel->SetMinEnergy(fMin);
  fModel->SetMaxEnergy(fMax);
  aP->AddDataSet(fAntiNucleonData);
  aP->RegisterMe(fModel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UrQMDAntiBarionBuilder::Build(G4AntiDeuteronInelasticProcess * aP)
{
  fModel->SetMinEnergy(fMin);
  fModel->SetMaxEnergy(fMax);
  aP->AddDataSet(fAntiNucleonData);
  aP->RegisterMe(fModel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UrQMDAntiBarionBuilder::Build(G4AntiTritonInelasticProcess * aP)
{
  fModel->SetMinEnergy(fMin);
  fModel->SetMaxEnergy(fMax);
  aP->AddDataSet(fAntiNucleonData);
  aP->RegisterMe(fModel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UrQMDAntiBarionBuilder::Build(G4AntiHe3InelasticProcess * aP)
{
  fModel->SetMinEnergy(fMin);
  fModel->SetMaxEnergy(fMax);
  aP->AddDataSet(fAntiNucleonData);
  aP->RegisterMe(fModel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UrQMDAntiBarionBuilder::Build(G4AntiAlphaInelasticProcess * aP)
{
  fModel->SetMinEnergy(fMin);
  fModel->SetMaxEnergy(fMax);
  aP->AddDataSet(fAntiNucleonData);
  aP->RegisterMe(fModel);
}
#endif //G4_USE_URQMD

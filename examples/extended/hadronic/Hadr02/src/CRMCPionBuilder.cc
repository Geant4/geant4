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
/// \file hadronic/Hadr02/src/CRMCPionBuilder.cc
/// \brief Implementation of the CRMCPionBuilder class
//
//
//---------------------------------------------------------------------------
//
// ClassName: CRMCPionBuilder
//
// Author:    2018 Alberto Ribon
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifdef G4_USE_CRMC

#include "CRMCPionBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4HadronicParameters.hh"
#include "G4SystemOfUnits.hh"


CRMCPionBuilder::CRMCPionBuilder() {
  fMin = 0.0*MeV;  // This value does not matter in practice because we are going
                   // to use this model only at high energies.
  fMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  fModel = new G4CRMCModel;
}


CRMCPionBuilder::~CRMCPionBuilder() {}


void CRMCPionBuilder::Build( G4HadronElasticProcess* ) {}


void CRMCPionBuilder::Build( G4PionPlusInelasticProcess* aP ) {
  fModel->SetMinEnergy( fMin );
  fModel->SetMaxEnergy( fMax );
  aP->RegisterMe( fModel );
}


void CRMCPionBuilder::Build( G4PionMinusInelasticProcess* aP ) {
  fModel->SetMinEnergy( fMin );
  fModel->SetMaxEnergy( fMax );
  aP->RegisterMe( fModel );
}

#endif //G4_USE_CRMC


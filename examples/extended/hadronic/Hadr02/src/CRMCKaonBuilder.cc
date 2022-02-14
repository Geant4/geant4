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
/// \file hadronic/Hadr02/src/CRMCKaonBuilder.cc
/// \brief Implementation of the CRMCKaonBuilder class
//
//
//---------------------------------------------------------------------------
//
// ClassName: CRMCKaonBuilder
//
// Author:    2018 Alberto Ribon
//
// Modified:
// -  18-May-2021 Alberto Ribon : Used the latest Geant4-CRMC interface.
//
//----------------------------------------------------------------------------
//
#ifdef G4_USE_CRMC

#include "CRMCKaonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4HadronInelasticProcess.hh"
#include "HadronicInelasticModelCRMC.hh"
#include "G4HadronicParameters.hh"
#include "G4SystemOfUnits.hh"


CRMCKaonBuilder::CRMCKaonBuilder( const G4int crmcModelId, const std::string & crmcModelName ) {
  fMin = 0.0*MeV;  // This value does not matter in practice because we are going
                   // to use this model only at high energies.
  fMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  fModel = new HadronicInelasticModelCRMC( crmcModelId, crmcModelName );
}


CRMCKaonBuilder::~CRMCKaonBuilder() {}


void CRMCKaonBuilder::Build( G4HadronElasticProcess* ) {}


void CRMCKaonBuilder::Build( G4HadronInelasticProcess* aP ) {
  fModel->SetMinEnergy( fMin );
  fModel->SetMaxEnergy( fMax );
  aP->RegisterMe( fModel );
}

#endif //G4_USE_CRMC

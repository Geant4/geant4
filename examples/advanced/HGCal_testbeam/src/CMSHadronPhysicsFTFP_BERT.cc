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
#include "CMSHadronPhysicsFTFP_BERT.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"

CMSHadronPhysicsFTFP_BERT::CMSHadronPhysicsFTFP_BERT(G4int)
    : G4HadronPhysicsFTFP_BERT("hInelastic FTFP_BERT", false) {
  minFTFP_pion = 3 * CLHEP::GeV; 
  maxBERT_pion = 12 * CLHEP::GeV;
  minFTFP_kaon = 3 * CLHEP::GeV;
  maxBERT_kaon = 6 * CLHEP::GeV;
  minFTFP_proton = 3 * CLHEP::GeV;
  maxBERT_proton = 6 * CLHEP::GeV;
  minFTFP_neutron = 3 * CLHEP::GeV;
  maxBERT_neutron = 6 * CLHEP::GeV;
}

CMSHadronPhysicsFTFP_BERT::~CMSHadronPhysicsFTFP_BERT() {}

void CMSHadronPhysicsFTFP_BERT::ConstructProcess() {
  if (G4Threading::IsMasterThread()) {
    DumpBanner();
  }
  CreateModels();
}

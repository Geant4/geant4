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
// $Id: HadrontherapyMuonStandard.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "HRMuonMinusCapture.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4MuonMinusCaptureAtRest.hh" 


HRMuonMinusCapture::HRMuonMinusCapture(const G4String& name): 
   G4VPhysicsConstructor(name)
{ 
  G4cout<< "HADRONIC PROCESS(ES): G4MuonMinusCaptureAtRest (muon-)" 
        << G4endl
        << "APPLIED MODEL(S): -" 
        << G4endl;
}

HRMuonMinusCapture::~HRMuonMinusCapture()
{ }

void HRMuonMinusCapture::ConstructProcess()
{

  // *************
  // *** Muon- ***
  // *************

  G4MuonMinusCaptureAtRest* muonMinusCaptureProcess= new G4MuonMinusCaptureAtRest();

  G4ParticleDefinition* particle = G4MuonMinus::MuonMinus(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddProcess(muonMinusCaptureProcess, 0, -1, -1);

}

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
// $Id: G4ThermalNeutrons.cc 70995 2013-06-09 00:56:34Z adotti $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4ThermalNeutrons
//
// Author: 4 May 2017 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4ThermalNeutrons.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4HadronicProcess.hh"

#include "G4ParticleHPThermalScattering.hh"
#include "G4ParticleHPThermalScatteringData.hh"

#include "G4BuilderType.hh"

#include "G4SystemOfUnits.hh"

G4ThermalNeutrons::G4ThermalNeutrons(G4int ver) :
  G4VHadronPhysics("G4ThermalNeutrons"), verbose(ver) {
}

G4ThermalNeutrons::~G4ThermalNeutrons() {}

void G4ThermalNeutrons::ConstructProcess() {

  if(verbose > 0) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4Neutron* part = G4Neutron::Neutron();
  G4HadronicProcess* hpel = FindElasticProcess(part);
  if(!hpel) {
    G4cout << "### " << GetPhysicsName() 
	   << " WARNING: Fail to add thermal neutron scattering" << G4endl;
    return;
  }

  G4int ni = (hpel->GetHadronicInteractionList()).size();
  if(ni < 1) {
    G4cout << "### " << GetPhysicsName() 
	   << " WARNING: Fail to add thermal neutron scattering - Nint= " 
	   << ni << G4endl;
    return;
  }
  (hpel->GetHadronicInteractionList())[ni-1]->SetMinEnergy(4*CLHEP::eV);

  hpel->RegisterMe(new G4ParticleHPThermalScattering());
  hpel->AddDataSet(new G4ParticleHPThermalScatteringData());
  
}

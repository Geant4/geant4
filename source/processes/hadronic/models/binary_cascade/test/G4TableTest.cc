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
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4TableTest.cc
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it), 
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"

#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <assert.h>

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

#include "Randomize.hh"

#include "G4ParticleDefinition.hh"
#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4Gamma.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonPlus.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4Lambda.hh"

#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"

#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"

#include "G4VShortLivedParticle.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedTable.hh"

#include "G4VXResonanceTable.hh"
#include "G4XNNstarTable.hh"
#include "G4XDeltaNstarTable.hh"
#include "G4XDeltaDeltaTable.hh"
#include "G4XNDeltaTable.hh"
#include "G4XNDeltastarTable.hh"
#include "G4XDeltaDeltastarTable.hh"

#include "G4ResonanceNames.hh"

#include <vector>

int main()
{
  // MGP ---- HBOOK initialization
  HepTupleManager* hbookManager;
  G4String fileName = "xtables.hbook";
  hbookManager = new HBookFile(fileName, 58);
  assert (hbookManager != 0);

  // ==== Initialization phase ====

  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();

  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
  G4ParticleDefinition* antiProton = G4AntiProton::AntiProtonDefinition();
  G4ParticleDefinition* neutron = G4Neutron::NeutronDefinition();   
  G4ParticleDefinition* antiNeutron = G4AntiNeutron::AntiNeutronDefinition();   

  G4ParticleDefinition* pionPlus = G4PionPlus::PionPlusDefinition();
  G4ParticleDefinition* pionMinus = G4PionMinus::PionMinusDefinition();
  G4ParticleDefinition* pionZero = G4PionZero::PionZeroDefinition();

  G4ParticleDefinition* kaonPlus = G4KaonPlus::KaonPlusDefinition();
  G4ParticleDefinition* kaonMinus = G4KaonMinus::KaonMinusDefinition();
   
  G4ParticleDefinition* lambda = G4Lambda::LambdaDefinition();

  G4ParticleDefinition* theMuonPlus = G4MuonPlus::MuonPlusDefinition();
  G4ParticleDefinition* theMuonMinus = G4MuonMinus::MuonMinusDefinition();

  G4ParticleDefinition* theNeutrinoMu = G4NeutrinoMu::NeutrinoMuDefinition();
  G4ParticleDefinition* theAntiNeutrinoMu = G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  // Construct resonances
  G4ShortLivedConstructor shortLived;
  shortLived.ConstructParticle();

  // Get the particle table
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  // ==== End of the initialization phase ====

  // The list of resonances handled by the Kinetic Model
  G4ResonanceNames* resonanceList = new G4ResonanceNames;

  // ==== Nucleons ====

  G4cout << proton->GetParticleName() 
	 <<": mass " << proton->GetPDGMass()
	 <<", width " << proton->GetPDGWidth()
	 <<", iSpin "<< proton->GetPDGiSpin()
	 <<", iIsospin " << proton->GetPDGiIsospin()
	 <<", iIsospin3 " << proton->GetPDGiIsospin3() 
	 << G4endl; 

  G4cout << neutron->GetParticleName() 
	 <<": mass " << neutron->GetPDGMass()
	 <<", width " << neutron->GetPDGWidth()
	 <<", iSpin "<< neutron->GetPDGiSpin()
	 <<", iIsospin " << neutron->GetPDGiIsospin()
	 <<", iIsospin3 " << neutron->GetPDGiIsospin3() 
	 << G4endl; 

  // ===== Delta =====

  std::vector<G4String> listDelta = resonanceList->DeltaNames();
  G4int nDelta = listDelta.size();
  G4cout << G4endl << "===== Delta ===== " << G4endl;

  G4int i;
  for (i=0; i<nDelta; i++)
    {
      // Particle information
      G4String name = listDelta[i];
      G4ParticleDefinition* def = particleTable->FindParticle(name);
      if (def == 0) G4cout << name << "does not have a ParticleDefinition " << G4endl;
      G4cout << def->GetParticleName() 
	     <<": mass " << def->GetPDGMass()
	     <<", width " << def->GetPDGWidth()
	     <<", iSpin "<< def->GetPDGiSpin()
	     <<", iIsospin " << def->GetPDGiIsospin()
	     <<", iIsospin3 " << def->GetPDGiIsospin3() 
	     << G4endl; 
    }

  // pp -> N Delta cross section table
  G4VXResonanceTable* tableND = new G4XNDeltaTable;
  
  for (i=0; i<nDelta; i+=4)
    {
      // Book a ntuple
      G4String name = listDelta[i];
      G4String ntName = "N " + name;
      HepTuple* ntupleND = hbookManager->ntuple(ntName);
      
      G4String xName = " p p -> Delta " + name;
      const G4PhysicsVector* sigma = tableND->CrossSectionTable(name);
      G4int entries = ((G4PhysicsVector*) sigma)->GetVectorLength();
      G4cout << "----------- " << xName <<" cross section table ----------- " 
	     << G4endl;

      // Fill the ntuple
      G4int j;
      for (j=0; j<entries; j++)
	{
	  const G4double energy = ((G4PhysicsVector*) sigma)->GetLowEdgeEnergy(j);
	  G4bool dummy = false;
	  G4double x = ((G4PhysicsVector*) sigma)->GetValue(energy,dummy) / millibarn;
	  G4cout << j << ") energy = " << energy /GeV << " GeV - sigma = " << x << " mb " 
		 << G4endl;
	  ntupleND->column("e",energy /GeV);
	  ntupleND->column("sigma",x);
	  ntupleND->dumpData();
	}
      G4cout << ntName << " ntuple available" << G4endl;
    }
  delete tableND;

  // pp -> Delta Delta cross section table
  G4VXResonanceTable* tableDD = new G4XDeltaDeltaTable;

  for (i=0; i<nDelta; i+=4)
    {
      // Book a ntuple
      G4String name = listDelta[i];
      G4String ntName = "delta " + name;
      HepTuple* ntupleDD = hbookManager->ntuple(ntName);
      
      G4String xName = " p p -> Delta " + name;
      const G4PhysicsVector* sigma = tableDD->CrossSectionTable(name);
      G4int entries = ((G4PhysicsVector*) sigma)->GetVectorLength();
      G4cout << "----------- " << xName << " cross section table ----------- " 
	     << G4endl;

      // Fill the ntuple
      G4int j;
      for (j=0; j<entries; j++)
	{
	  const G4double energy = ((G4PhysicsVector*) sigma)->GetLowEdgeEnergy(j);
	  G4bool dummy = false;
	  G4double x = ((G4PhysicsVector*) sigma)->GetValue(energy,dummy) / millibarn;
	  G4cout << j << ") energy = " << energy /GeV << " GeV - sigma = " << x << " mb " 
		 << G4endl;
	  ntupleDD->column("e",energy /GeV);
	  ntupleDD->column("sigma",x);
	  ntupleDD->dumpData();
	}
      G4cout << ntName << " ntuple available" << G4endl;
    }
  delete tableDD;


  // Excited Nucleons =====

  G4cout << G4endl << "===== Excited Nucleons ===== " << G4endl;

  std::vector<G4String> listNstar = resonanceList->NstarNames();
  G4int nNstar = listNstar.size();

  for (i=0; i<nNstar; i++)
    {
      G4String name = listNstar[i];
      G4ParticleDefinition* def = particleTable->FindParticle(name);
      if (def == 0) G4cout << name << "does not have a ParticleDefinition " << G4endl;
      G4cout << def->GetParticleName() 
	     <<": mass " << def->GetPDGMass()
	     <<", width " << def->GetPDGWidth()
	     <<", iSpin "<< def->GetPDGiSpin()
	     <<", iIsospin " << def->GetPDGiIsospin()
	     <<", iIsospin3 " << def->GetPDGiIsospin3() 
	     << G4endl; 
    }

  // pp -> N Nstar  cross section table
  G4VXResonanceTable* tableNNstar = new G4XNNstarTable;
  
  for (i=0; i<nNstar; i+=2)
    {
      // Book a ntuple
      G4String name = listNstar[i];
      G4String ntName = "N " + name;
      HepTuple* ntupleNNstar = hbookManager->ntuple(ntName);

      G4String xName = " p p -> N " + name;
      const G4PhysicsVector* sigma = tableNNstar->CrossSectionTable(name);
      if (sigma == 0)
	G4cout << name << " does not return a valid PhysicsVector " << G4endl;

      G4int entries = ((G4PhysicsVector*) sigma)->GetVectorLength();
      G4cout << "----------- " << xName << " cross section table ----------- " 
	     << G4endl;

      // Fill the ntuple
      G4int j;
      for (j=0; j<entries; j++)
	{
	  const G4double energy = ((G4PhysicsVector*) sigma)->GetLowEdgeEnergy(j);
	  G4bool dummy = false;
	  G4double x = ((G4PhysicsVector*) sigma)->GetValue(energy,dummy) / millibarn;
	  G4cout << j << ") energy = " << energy /GeV << " GeV - sigma = " << x << " mb " 
		 << G4endl;
	  ntupleNNstar->column("e",energy /GeV);
	  ntupleNNstar->column("sigma",x);
	  ntupleNNstar->dumpData();
	}
      G4cout << ntName << " ntuple available" << G4endl;
    }
  delete tableNNstar;

  // pp -> Nstar Delta cross section table
  G4VXResonanceTable* tableDNstar = new G4XDeltaNstarTable;

  for (i=0; i<nNstar; i+=2)
    {
      // Book a ntuple
      G4String name = listNstar[i];
      G4String ntName = "delta " + name;
      HepTuple* ntupleDNstar = hbookManager->ntuple(ntName);

      G4String xName = " p p -> Delta " + name;
      const G4PhysicsVector* sigma = tableDNstar->CrossSectionTable(name);
      G4int entries = ((G4PhysicsVector*) sigma)->GetVectorLength();
      G4cout << "----------- " << xName << " cross section table ----------- " 
	     << G4endl;

      // Fill the ntuple
      G4int j;
      for (j=0; j<entries; j++)
	{
	  const G4double energy = ((G4PhysicsVector*) sigma)->GetLowEdgeEnergy(j);
	  G4bool dummy = false;
	  G4double x = ((G4PhysicsVector*) sigma)->GetValue(energy,dummy) / millibarn;
	  G4cout << j << ") energy = " << energy /GeV << " GeV - sigma = " << x << " mb " 
		 << G4endl;
	  ntupleDNstar->column("e",energy /GeV);
	  ntupleDNstar->column("sigma",x);
	  ntupleDNstar->dumpData();
	}
      G4cout << ntName << " ntuple available" << G4endl;
    }
  delete tableDNstar;

  // ===== Excited Deltas ===== 

  G4cout << G4endl << "===== Excited Deltas ===== " << G4endl;

  std::vector<G4String> listDeltastar = resonanceList->DeltastarNames();
  G4int nDeltastar = listDeltastar.size();

  for (i=0; i<nDeltastar; i++)
    {
      G4String name = listDeltastar[i];
      G4ParticleDefinition* def = particleTable->FindParticle(name);
      if (def == 0) G4cout << name << "does not have a ParticleDefinition " << G4endl;
      G4cout << def->GetParticleName() 
	     <<": mass " << def->GetPDGMass()
	     <<", width " << def->GetPDGWidth()
	     <<", iSpin "<< def->GetPDGiSpin()
	     <<", iIsospin " << def->GetPDGiIsospin()
	     <<", iIsospin3 " << def->GetPDGiIsospin3() 
	     << G4endl; 
    }
  
  // pp -> N Deltastar cross section table
  G4VXResonanceTable* tableNDstar = new G4XNDeltastarTable;
       
  for (i=0; i<nDeltastar; i+=4)
    {
      // Book a ntuple
      G4String name = listDeltastar[i];
      G4String ntName = "N " + name;
      HepTuple* ntupleNDstar = hbookManager->ntuple(ntName);

      G4String xName = " p p -> N " + name;
      const G4PhysicsVector* sigma = tableNDstar->CrossSectionTable(name);
      G4int entries = ((G4PhysicsVector*) sigma)->GetVectorLength();
      G4cout << "----------- " << xName << " cross section table ----------- " 
	     << G4endl;

  // Fill the ntuple
      G4int j;
      for (j=0; j<entries; j++)
	{
	  const G4double energy = ((G4PhysicsVector*) sigma)->GetLowEdgeEnergy(j);
	  G4bool dummy = false;
	  G4double x = ((G4PhysicsVector*) sigma)->GetValue(energy,dummy) / millibarn;
	  G4cout << j << ") energy = " << energy /GeV << " GeV - sigma = " << x << " mb " 
		 << G4endl;
	  ntupleNDstar->column("e",energy /GeV);
	  ntupleNDstar->column("sigma",x);
	  ntupleNDstar->dumpData();
	}
      G4cout << ntName << " ntuple available" << G4endl;
    }
  delete tableNDstar;

  // pp -> Delta Deltastar cross section table
  G4VXResonanceTable* tableDDstar = new G4XDeltaDeltastarTable;
      
  for (i=0; i<nDeltastar; i+=4)
    {
      // Book a ntuple
      G4String name = listDeltastar[i];
      G4String ntName = "delta " + name;
      HepTuple* ntupleDDstar = hbookManager->ntuple(ntName);

      G4String xName = " p p -> Delta " + name;
      const G4PhysicsVector* sigma = tableDDstar->CrossSectionTable(name);
      G4int entries = ((G4PhysicsVector*) sigma)->GetVectorLength();
      G4cout << "----------- " << xName << " cross section table ----------- " 
	     << G4endl;

      // Fill the ntuple
      G4int j;
      for (j=0; j<entries; j++)
	{
	  const G4double energy = ((G4PhysicsVector*) sigma)->GetLowEdgeEnergy(j);
	  G4bool dummy = false;
	  G4double x = ((G4PhysicsVector*) sigma)->GetValue(energy,dummy) / millibarn;
	  G4cout << j << ") energy = " << energy /GeV << " GeV - sigma = " << x << " mb " 
		 << G4endl;
	  ntupleDDstar->column("e",energy /GeV);
	  ntupleDDstar->column("sigma",x);
	  ntupleDDstar->dumpData();
	}
      G4cout << ntName << " ntuple available" << G4endl;
    }
  delete tableDDstar;

  hbookManager->write();

  G4cout << "ntuples are in file " << fileName<< G4endl;

  delete resonanceList;

  return EXIT_SUCCESS;
}

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
/// \file electromagnetic/TestEm0/DirectAccess.cc
/// \brief Main program of the electromagnetic/TestEm0 example
//
//
// 
// ------------------------------------------------------------
//
//  To print cross sections per atom and mean free path for simple material
//
#include "G4Material.hh"

#include "G4PEEffectFluoModel.hh"
#include "G4KleinNishinaCompton.hh"
#include "G4BetheHeitlerModel.hh"

#include "G4eeToTwoGammaModel.hh"

#include "G4MollerBhabhaModel.hh"
#include "G4SeltzerBergerModel.hh"

#include "G4BetheBlochModel.hh"
#include "G4BraggModel.hh"

#include "G4MuBetheBlochModel.hh"
#include "G4MuBremsstrahlungModel.hh"
#include "G4MuPairProductionModel.hh"

#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4MuonPlus.hh"

#include "G4DataVector.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"

int main() {

  G4UnitDefinition::BuildUnitsTable();

  G4ParticleDefinition* gamma = G4Gamma::Gamma();
  G4ParticleDefinition* posit = G4Positron::Positron();
  G4ParticleDefinition* elec = G4Electron::Electron();
  G4ParticleDefinition* prot = G4Proton::Proton();
  G4ParticleDefinition* muon = G4MuonPlus::MuonPlus();
  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  G4DataVector cuts;
  cuts.push_back(1*keV);

  // define materials
  //
  G4Material* material = 
    G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  G4MaterialCutsCouple* couple = new G4MaterialCutsCouple(material);
  couple->SetIndex(0);

  // work only for simple materials
  G4double Z = material->GetZ();
  G4double A = material->GetA();

  // initialise gamma processes (models)
  //  
  G4VEmModel* phot = new G4PEEffectFluoModel();
  G4VEmModel* comp = new G4KleinNishinaCompton();
  G4VEmModel* conv = new G4BetheHeitlerModel(); 
  phot->Initialise(gamma, cuts);
  comp->Initialise(gamma, cuts);
  conv->Initialise(gamma, cuts);

  // valid pointer to a couple is needed for this model
  phot->SetCurrentCouple(couple);

  // compute CrossSection per atom and MeanFreePath
  //
  G4double Emin = 1.01*MeV, Emax = 2.01*MeV, dE = 100*keV;

  G4cout << "\n #### Gamma : CrossSectionPerAtom and MeanFreePath for " 
         << material->GetName() << G4endl;
  G4cout << "\n Energy \t PhotoElec \t Compton \t Conversion \t";
  G4cout <<           "\t PhotoElec \t Compton \t Conversion" << G4endl;
  
  for (G4double Energy = Emin; Energy <= Emax; Energy += dE) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
     << "\t" 
     << G4BestUnit (phot->ComputeCrossSectionPerAtom(gamma,Energy,Z),"Surface")
     << "\t"         
     << G4BestUnit (comp->ComputeCrossSectionPerAtom(gamma,Energy,Z),"Surface")
     << "\t"         
     << G4BestUnit (conv->ComputeCrossSectionPerAtom(gamma,Energy,Z),"Surface")
     << "\t \t"         
     << G4BestUnit (phot->ComputeMeanFreePath(gamma,Energy,material),"Length")
     << "\t"         
     << G4BestUnit (comp->ComputeMeanFreePath(gamma,Energy,material),"Length")
     << "\t"         
     << G4BestUnit (conv->ComputeMeanFreePath(gamma,Energy,material),"Length");
  }

  G4cout << G4endl;

  // initialise positron annihilation (model)
  //    
  G4VEmModel* anni = new G4eeToTwoGammaModel();
  anni->Initialise(posit, cuts);
  
  // compute CrossSection per atom and MeanFreePath
  //
  Emin = 1.01*MeV; Emax = 2.01*MeV; dE = 100*keV;

  G4cout << "\n #### e+ annihilation : CrossSectionPerAtom and MeanFreePath"
         << " for " << material->GetName() << G4endl;
  G4cout << "\n Energy \t e+ annihil \t";
  G4cout <<           "\t e+ annihil" << G4endl;
  
  for (G4double Energy = Emin; Energy <= Emax; Energy += dE) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
     << "\t" 
     << G4BestUnit (anni->ComputeCrossSectionPerAtom(posit,Energy,Z),"Surface")
     << "\t \t"         
     << G4BestUnit (anni->ComputeMeanFreePath(posit,Energy,material),"Length");
  }

  G4cout << G4endl;

  // initialise electron processes (models)
  //    
  G4VEmModel* ioni = new G4MollerBhabhaModel();
  G4VEmModel* brem = new G4SeltzerBergerModel();
  ioni->Initialise(elec, cuts);
  brem->Initialise(elec, cuts);

  // compute CrossSection per atom and MeanFreePath
  //
  Emin = 1.01*MeV; Emax = 101.01*MeV; dE = 10*MeV;
  G4double Ecut = 100*keV;

  G4cout << "\n ####electron: CrossSection, MeanFreePath and StoppingPower"
         << " for " << material->GetName() 
         << ";\tEnergy cut = " << G4BestUnit (Ecut, "Energy") << G4endl;
         
  G4cout << "\n Energy \t ionization \t bremsstra \t";
  G4cout <<           "\t ionization \t bremsstra \t";
  G4cout <<           "\t ionization \t bremsstra" << G4endl;
  
  for (G4double Energy = Emin; Energy <= Emax; Energy += dE) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
     << "\t" 
     << G4BestUnit (ioni->ComputeCrossSectionPerAtom(elec,Energy,Z,A,Ecut),
                   "Surface")
     << "\t" 
     << G4BestUnit (brem->ComputeCrossSectionPerAtom(elec,Energy,Z,A,Ecut),
                   "Surface")                   
     << "\t \t"         
     << G4BestUnit (ioni->ComputeMeanFreePath(elec,Energy,material,Ecut),
                   "Length")
     << "\t"         
     << G4BestUnit (brem->ComputeMeanFreePath(elec,Energy,material,Ecut),
                   "Length")                   
     << "\t \t"         
     << G4BestUnit (ioni->ComputeDEDXPerVolume(material,elec,Energy,Ecut),
                   "Energy/Length")
     << "\t"         
     << G4BestUnit (brem->ComputeDEDXPerVolume(material,elec,Energy,Ecut),
                   "Energy/Length");                                      
  }
  
  G4cout << G4endl;

  // initialise proton processes (models)
  //    
  ioni = new G4BetheBlochModel();
  ioni->Initialise(prot, cuts);
  
  // compute CrossSection per atom and MeanFreePath
  //
  Emin = 1.01*MeV; Emax = 102.01*MeV; dE = 10*MeV;
  Ecut = 100*keV;

  G4cout << "\n #### proton : CrossSection, MeanFreePath and StoppingPower"
         << " for " << material->GetName() 
         << ";\tEnergy cut = " << G4BestUnit (Ecut, "Energy") << G4endl;
         
  G4cout << "\n Energy \t ionization \t";
  G4cout <<           "\t ionization \t";
  G4cout <<           "\t ionization" << G4endl;
  
  for (G4double Energy = Emin; Energy <= Emax; Energy += dE) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
     << "\t" 
     << G4BestUnit (ioni->ComputeCrossSectionPerAtom(prot,Energy,Z,A,Ecut),
                   "Surface")
     << "\t \t"         
     << G4BestUnit (ioni->ComputeMeanFreePath(prot,Energy,material,Ecut),
                   "Length")           
     << "\t \t"         
     << G4BestUnit (ioni->ComputeDEDXPerVolume(material,prot,Energy,Ecut),
                   "Energy/Length");
  }
  
  G4cout << G4endl;
  
  // low energy : Bragg Model
  ioni = new G4BraggModel(prot);
  ioni->Initialise(prot, cuts);
  
  // compute CrossSection per atom and MeanFreePath
  //
  Emin = 1.1*keV; Emax = 2.01*MeV; dE = 300*keV;
  Ecut = 10*keV;
  
  G4cout << "\n #### proton : low energy model (Bragg) "
         << ";\tEnergy cut = " << G4BestUnit (Ecut, "Energy") << G4endl;
                           
  G4cout << "\n Energy \t ionization \t";
  G4cout <<           "\t ionization \t";
  G4cout <<           "\t ionization" << G4endl;
  
  for (G4double Energy = Emin; Energy <= Emax; Energy += dE) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
     << "\t" 
     << G4BestUnit (ioni->ComputeCrossSectionPerAtom(prot,Energy,Z,A,Ecut),
                   "Surface")
     << "\t \t"         
     << G4BestUnit (ioni->ComputeMeanFreePath(prot,Energy,material,Ecut),
                   "Length")           
     << "\t \t"         
     << G4BestUnit (ioni->ComputeDEDXPerVolume(material,prot,Energy,Ecut),
                   "Energy/Length");
  }
  
  G4cout << G4endl;
  
  // initialise muon processes (models)
  //  
  ioni = new G4MuBetheBlochModel();
  brem = new G4MuBremsstrahlungModel();
  G4VEmModel* pair = new G4MuPairProductionModel();
  ioni->Initialise(muon, cuts);
  brem->Initialise(muon, cuts);
  pair->Initialise(muon, cuts);
   
  // compute CrossSection per atom and MeanFreePath
  //
  Emin = 1.01*GeV; Emax = 101.01*GeV; dE = 10*GeV;
  Ecut = 10*MeV;

  G4cout << "\n ####muon: CrossSection and MeanFreePath for "
         << material->GetName() 
         << ";\tEnergy cut = " << G4BestUnit (Ecut, "Energy") << G4endl;
         
  G4cout << "\n Energy \t ionization \t bremsstra \t pair_prod \t";
  G4cout <<           "\t ionization \t bremsstra \t pair_prod" << G4endl;
  
  for (G4double Energy = Emin; Energy <= Emax; Energy += dE) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
     << "\t" 
     << G4BestUnit (ioni->ComputeCrossSectionPerAtom(muon,Energy,Z,A,Ecut),
                   "Surface")
     << "\t" 
     << G4BestUnit (brem->ComputeCrossSectionPerAtom(muon,Energy,Z,A,Ecut),
                   "Surface")
     << "\t"                    
     << G4BestUnit (pair->ComputeCrossSectionPerAtom(muon,Energy,Z,A,Ecut),
                   "Surface")                                      
     << "\t \t"         
     << G4BestUnit (ioni->ComputeMeanFreePath(muon,Energy,material,Ecut),
                   "Length")
     << "\t"         
     << G4BestUnit (brem->ComputeMeanFreePath(muon,Energy,material,Ecut),
                   "Length")
     << "\t"         
     << G4BestUnit (pair->ComputeMeanFreePath(muon,Energy,material,Ecut),
                   "Length");
  }
  
  G4cout << G4endl;
  
  G4cout << "\n ####muon: StoppingPower for "
         << material->GetName() 
         << ";\tEnergy cut = " << G4BestUnit (Ecut, "Energy") << G4endl;
         
  G4cout << "\n Energy \t ionization \t bremsstra \t pair_prod \t" << G4endl;
  
  for (G4double Energy = Emin; Energy <= Emax; Energy += dE) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
     << "\t"         
     << G4BestUnit (ioni->ComputeDEDXPerVolume(material,muon,Energy,Ecut),
                   "Energy/Length")
     << "\t"         
     << G4BestUnit (brem->ComputeDEDXPerVolume(material,muon,Energy,Ecut),
                   "Energy/Length")
     << "\t"         
     << G4BestUnit (pair->ComputeDEDXPerVolume(material,muon,Energy,Ecut),
                   "Energy/Length");
  }
  
  G4cout << G4endl;                                   
  return EXIT_SUCCESS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

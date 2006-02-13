//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: CrossSectionV52.cc,v 1.3 2006-02-13 16:52:24 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// ------------------------------------------------------------
//
//  To print cross sections per atom and mean free path for simple material
//
#include "G4Material.hh"

#include "G4PhotoElectricEffect52.hh"
#include "G4ComptonScattering52.hh"
#include "G4GammaConversion52.hh"

#include "G4eplusAnnihilation52.hh"

#include "G4eIonisation52.hh"
#include "G4eBremsstrahlung52.hh"

#include "G4hIonisation52.hh"

#include "G4MuIonisation52.hh"
#include "G4MuBremsstrahlung52.hh"
#include "G4MuPairProduction52.hh"

#include "globals.hh"
#include "G4UnitsTable.hh"

int main() {

  G4UnitDefinition::BuildUnitsTable();

  // define materials
  //
  G4double Z, A;

  G4Material* material =
  new G4Material("Iodine", Z=53., A=126.90*g/mole, 4.93*g/cm3);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
   
  // initialise gamma processes
  // 
  G4PhotoElectricEffect52* phot = new G4PhotoElectricEffect52();
  G4ComptonScattering52*   comp = new G4ComptonScattering52();
  G4GammaConversion52*     conv = new G4GammaConversion52();
  
  // compute CrossSection per atom and MeanFreePath
  //
  G4double Emin = 1.01*MeV, Emax = 2.01*MeV, dE = 100*keV;

  G4cout << "\n #### Gamma CrossSectionPerAtom and MeanFreePath for " 
         << material->GetName() << G4endl;
  G4cout << "\n Energy \t PhotoElec \t Compton \t Conversion \t";
  G4cout <<           "\t PhotoElec \t Compton \t Conversion" << G4endl;
  
  for (G4double Energy = Emin; Energy <= Emax; Energy += dE) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
         << "\t" 
         << G4BestUnit (phot->ComputeCrossSectionPerAtom(Energy,Z), "Surface")
         << "\t"	 
	 << G4BestUnit (comp->ComputeCrossSectionPerAtom(Energy,Z), "Surface")
         << "\t"	 
	 << G4BestUnit (conv->ComputeCrossSectionPerAtom(Energy,Z), "Surface")
         << "\t \t"	 
	 << G4BestUnit (phot->ComputeMeanFreePath(Energy,material), "Length")
         << "\t"	 
	 << G4BestUnit (comp->ComputeMeanFreePath(Energy,material), "Length")
         << "\t"	 
	 << G4BestUnit (conv->ComputeMeanFreePath(Energy,material), "Length");	 
  }

  G4cout << G4endl;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  // initialise positron annihilation
  // 
  G4eplusAnnihilation52* anni = new G4eplusAnnihilation52();
  
  // compute CrossSection per atom and MeanFreePath
  //
  Emin = 1.01*MeV; Emax = 2.01*MeV; dE = 100*keV;

  G4cout << "\n #### e+ annihilation CrossSectionPerAtom and MeanFreePath for " 
         << material->GetName() << G4endl;
  G4cout << "\n Energy \t e+ annihil \t";
  G4cout <<           "\t e+ annihil" << G4endl;
  
  for (G4double Energy = Emin; Energy <= Emax; Energy += dE) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
         << "\t" 
         << G4BestUnit (anni->ComputeCrossSectionPerAtom(Energy,Z), "Surface")
         << "\t \t"	 
	 << G4BestUnit (anni->ComputeMeanFreePath(Energy,material), "Length"); 
  }
  
  G4cout << G4endl;
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  // initialise electron processes
  // 
  G4ParticleDefinition* elec = G4Electron::Electron();
   
  G4eIonisation52* ioni = new G4eIonisation52();
  G4eBremsstrahlung52* brem = new G4eBremsstrahlung52();
  
  // compute CrossSection per atom and restricted dE/dx
  //
  Emin = 1.01*MeV; Emax = 101.01*MeV; dE = 10*MeV;
  G4double Ecut = 100*keV;

  G4cout << "\n ####electron: CrossSectionPerAtom and StoppingPower for "
         << material->GetName() 
	 << ";\tEnergy cut = " << G4BestUnit (Ecut, "Energy") << G4endl;
	 
  G4cout << "\n Energy \t ionization \t bremsstra \t";
  G4cout <<           "\t ionization" << G4endl;
  
  for (G4double Energy = Emin; Energy <= Emax; Energy += dE) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
     << "\t" 
     << G4BestUnit (ioni->ComputeCrossSectionPerAtom(*elec,Energy,Z,Ecut),
                   "Surface")
     << "\t" 
     << G4BestUnit (brem->ComputeCrossSectionPerAtom( elec,Energy,Z,Ecut),
                   "Surface")		   	   
     << "\t \t"	 
     << G4BestUnit (ioni->ComputeRestrictedMeandEdx(*elec,Energy,material,Ecut),
                   "Energy/Length");   		   
  }
  
  G4cout << G4endl;
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  // initialise hadron processes
  // 
  G4ParticleDefinition* prot = G4Proton::Proton();
   
  G4hIonisation52* ionis = new G4hIonisation52();
  
  // compute CrossSection per atom and restricted dE/dx
  //
  Emin = 1.01*MeV; Emax = 101.01*MeV; dE = 10*MeV;
  Ecut = 100*keV;

  G4cout << "\n ####proton: CrossSectionPerAtom and StoppingPower for "
         << material->GetName() 
	 << ";\tEnergy cut = " << G4BestUnit (Ecut, "Energy") << G4endl;
	 
  G4cout << "\n Energy \t ionization \t";
  G4cout <<           "\t ionization" << G4endl;
  
  for (G4double Energy = Emin; Energy <= Emax; Energy += dE) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
     << "\t" 
     << G4BestUnit (ionis->ComputeCrossSectionPerAtom(*prot,Energy,Z,Ecut),
                   "Surface")	   	   
     << "\t \t"	 
     << G4BestUnit (ionis->ComputeRestrictedMeandEdx(*prot,Energy,material,Ecut),
                   "Energy/Length");   		   
  }
  
  G4cout << G4endl;
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  // initialise muon processes
  // 
  G4ParticleDefinition* muon = G4MuonPlus::MuonPlus();
   
  G4MuIonisation52* muioni = new G4MuIonisation52();
  G4MuBremsstrahlung52* mubrem = new G4MuBremsstrahlung52();
  G4MuPairProduction52* mupair = new G4MuPairProduction52();
      
  // compute CrossSection per atom and restricted dE/dx
  //
  Emin = 1.01*GeV; Emax = 101.01*GeV; dE = 10*GeV;
  Ecut = 10*MeV;

  G4cout << "\n ####muon: CrossSectionPerAtom and StoppingPower for "
         << material->GetName() 
	 << ";\tEnergy cut = " << G4BestUnit (Ecut, "Energy") << G4endl;
	 
  G4cout << "\n Energy \t ionization \t bremsstra \t pair_prod  \t";
  G4cout <<           "\t ionization" << G4endl;
  
  for (G4double Energy = Emin; Energy <= Emax; Energy += dE) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
     << "\t" 
     << G4BestUnit (muioni->ComputeCrossSectionPerAtom(*muon,Energy,Z,Ecut),
                   "Surface")
     << "\t" 		   
     << G4BestUnit (mubrem->ComputeMicroscopicCrossSection(muon,Energy,Z,A,Ecut),
                   "Surface")
     << "\t" 		   
     << G4BestUnit (mupair->ComputeMicroscopicCrossSection(muon,Energy,Z,Ecut,
                                                                         Ecut),
                   "Surface")		   		   	   	   
     << "\t \t"	 
     << G4BestUnit (muioni->ComputeRestrictedMeandEdx(*muon,Energy,material,Ecut),
                   "Energy/Length");   		   
  }
  
  G4cout << G4endl;
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
                               
return EXIT_SUCCESS;
}

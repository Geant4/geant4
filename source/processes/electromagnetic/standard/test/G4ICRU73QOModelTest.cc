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
// 
// ------------------------------------------------------------
//
//  To print eloss by using of ICRU73QOModel
//
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4BetheBlochModel.hh"
#include "G4BraggModel.hh"
#include "G4ICRU73QOModel.hh"
#include "G4VEmModel.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleTable.hh"
#include "G4hIonisation.hh"
#include "G4PhysicsLogVector.hh"

#include "globals.hh"
#include "G4UnitsTable.hh"

#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

int main() {

  G4UnitDefinition::BuildUnitsTable();

  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial("G4_H");
  
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();

  // initialise proton processes (models)
  // 
  G4ParticleDefinition* part = G4AntiProton::AntiProton();
  //G4ParticleDefinition* part = G4MuonMinus::MuonMinus();

  G4cout << "pName = " << part->GetParticleName() << G4endl;
  G4String pName = part->GetParticleName();

  partTable->SetReadiness();

  G4VEmModel* bragg = new G4BraggModel();
  G4VEmModel* bethe = new G4BetheBlochModel();
  G4VEmModel* qo = new G4ICRU73QOModel();
 
  G4DataVector v;
  v.resize(100);

  bragg->Initialise(part,v);
  bethe->Initialise(part,v);
  qo->Initialise(part,v);
  
  G4double Emin = 1.1*keV; 
  G4double Emax = 100.01*MeV; 
 G4double Ecut = 1*GeV;
 
  size_t nBinTab = 35;

  G4PhysicsLogVector* pVector;
  pVector = new G4PhysicsLogVector(Emin, Emax, nBinTab);

  G4double Energy = 0.;
  G4double ELos1 = 0.;
  G4double ELos2 = 0.;
  G4double ELos3 = 0.;
   
  G4String mat = material->GetName();
  G4String asciiFileName = pName + "_QOEloss_diff_" + mat + ".out"; 

  std::ofstream asciiFile(asciiFileName, std::ios::out);
  if(asciiFile.is_open()) {
    asciiFile << " Energy(Mev)\t J(MeV/mm) (bragg  bethe  quantosc)" << G4endl;
  } else {
    G4cout << "ERROR file <" << asciiFileName << "> is not opened" << G4endl;
    exit(1);
  }

// Write to file
  
  for (size_t i = 0; i < nBinTab+1; i++) {
    Energy = pVector->GetLowEdgeEnergy(i); 
    ELos1   = bragg->ComputeDEDXPerVolume(material,part,Energy,Ecut);     
    ELos2   = bethe->ComputeDEDXPerVolume(material,part,Energy,Ecut);     
    ELos3   = qo->ComputeDEDXPerVolume(material,part,Energy,Ecut);     
    asciiFile << std::setiosflags(std::ios::fixed)
	      << std::setprecision(5)
	      << std::setiosflags(std::ios::right)
	      << std::setw(10);
    asciiFile << Energy;
    asciiFile << "\t";
    asciiFile << std::setiosflags(std::ios::fixed)
	      << std::setprecision(5)
	      << std::setiosflags(std::ios::right)
	      << std::setw(10);
    asciiFile << ELos1;
    asciiFile << "\t";
    asciiFile << std::setiosflags(std::ios::fixed)
	      << std::setprecision(5)
	      << std::setiosflags(std::ios::right)
	      << std::setw(10);
    asciiFile << ELos2;
    asciiFile << "\t";
    asciiFile << std::setiosflags(std::ios::fixed)
	      << std::setprecision(5)
	      << std::setiosflags(std::ios::right)
	      << std::setw(10);
    asciiFile << ELos3
	      << G4endl;
  }

    delete pVector;
  
    return EXIT_SUCCESS;
}

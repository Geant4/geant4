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

#include <vector>

int main() {
 
  G4UnitDefinition::BuildUnitsTable();

  G4NistManager* nist = G4NistManager::Instance();
  G4Material* material = nist->FindOrBuildMaterial("G4_H");
  
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  const std::vector<G4String> mnames = nist->GetNistMaterialNames();

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();

  // initialise proton processes (models)
  //   
  G4ParticleDefinition* part = G4AntiProton::AntiProton();
  G4ParticleDefinition* prot = G4Proton::Proton();
  //G4ParticleDefinition* part = G4MuonMinus::MuonMinus();

  G4cout << "pName = " << part->GetParticleName() << G4endl;
  G4String pName = part->GetParticleName();

  partTable->SetReadiness();

  G4VEmModel* bragg = new G4BraggModel();
  G4VEmModel* bethe1 = new G4BetheBlochModel();
  G4VEmModel* bethe2 = new G4BetheBlochModel();
  G4VEmModel* qo = new G4ICRU73QOModel();
 
  G4DataVector v;
  v.resize(100);

  bragg->Initialise(prot,v);
  bethe1->Initialise(prot,v);
  bethe2->Initialise(part,v);
  qo->Initialise(part,v);
  
  G4double Emin = 1.1*keV; 
  G4double Emax = 100.01*MeV; 
  G4double Ecut = 1*GeV;
 
  size_t nBinTab = 35;

  G4PhysicsLogVector* pVector;
  pVector = new G4PhysicsLogVector(Emin, Emax, nBinTab);

  G4double Energy = 0.;
  G4double ELos0 = 0.;
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
  G4double fact1[98];
  G4double fact2[98];
  G4double fact3[98]; 

  Energy = 2.*MeV;
  part = G4AntiProton::AntiProton();
  pName = part->GetParticleName(); 
  G4cout << "### E(MeV)= " << Energy/MeV << "  for " << pName 
	 << std::setprecision(4)
	 << G4endl;
  G4cout << " const G4double factorBethe[98] = { 1.0, " << G4endl;
  for(G4int j=0; j<98; ++j) {
    const G4Material* mat = nist->FindOrBuildMaterial(mnames[j]);
    G4int Z = G4int(mat->GetZ());
    ELos0   = bragg->ComputeDEDXPerVolume(mat,prot,Energy,Ecut);     
    ELos1   = bethe1->ComputeDEDXPerVolume(mat,prot,Energy,Ecut);     
    ELos2   = bethe2->ComputeDEDXPerVolume(mat,part,Energy,Ecut);     
    ELos3   = qo->ComputeDEDXPerVolume(mat,part,Energy,Ecut);
    fact1[j]= ELos2/ELos3;
    if(ELos2 > ELos0) {  fact1[j] *= ELos0/ELos1; }
    fact2[j]= ELos3*fact1[j]/ELos0;
    fact3[j]= ELos1/ELos2;
    G4cout << fact1[j] << ", ";
    if((Z/10)*10 == Z) { G4cout << "  // " << Z-9 << " - " << Z << G4endl; }
  } 
  G4cout << " } " << G4endl; 
  G4cout << " ##### ICRU73(pbar)/Bragg(p) ratio at 2 MeV: " << G4endl; 
  for(G4int j=1; j<98; ++j) {
    G4int Z = j;
    G4cout << fact2[j] << ", ";
    if((Z/10)*10 == Z) { G4cout << "  // " << Z-9 << " - " << Z << G4endl; }
  }
  G4cout << "  " << G4endl; 
  G4cout << " ##### p/pbar BetheBloch ratio at 2 MeV: " << G4endl; 
  for(G4int j=1; j<98; ++j) {
    G4int Z = j;
    G4cout << fact3[j] << ", ";
    if((Z/10)*10 == Z) { G4cout << "  // " << Z-9 << " - " << Z << G4endl; }
  }
  G4cout << "  " << G4endl; 
  part = G4Proton::Proton();
  pName = part->GetParticleName();
  G4cout << "### E(MeV)= " << Energy/MeV << "  for " << pName 
	 << " Bethe/Bragg" <<G4endl; 
  G4cout << " const G4double factorBragg[98] = { 1.0, " << G4endl;
  for(G4int j=0; j<98; ++j) {
    const G4Material* mat = nist->FindOrBuildMaterial(mnames[j]);
    G4int Z = G4int(mat->GetZ());
    ELos1   = bragg->ComputeDEDXPerVolume(mat,prot,Energy,Ecut);     
    ELos2   = bethe1->ComputeDEDXPerVolume(mat,prot,Energy,Ecut);     
    fact2[j]= ELos2/ELos1; 
    G4cout << fact2[j] << ", ";  
    if((Z/10)*10 == Z) { G4cout << "  // " << Z-9 << " - " << Z << G4endl; }
  }  
  G4cout << " } " << G4endl;  

  Energy = 2.*MeV;
  G4cout << "### E(MeV)= " << Energy/MeV << "  for " << pName 
	 << std::setprecision(5) << G4endl;
  for(G4int j=0; j<98; ++j) {
    const G4Material* mat = nist->FindOrBuildMaterial(mnames[j]);
    G4double Z = mat->GetZ();
    ELos1   = bragg->ComputeDEDXPerVolume(mat,prot,Energy,Ecut);     
    ELos2   = bethe1->ComputeDEDXPerVolume(mat,prot,Energy,Ecut);     
    ELos3   = qo->ComputeDEDXPerVolume(mat,part,Energy,Ecut);
    G4cout << "Z= " << Z << " p bragg/bethe-1= " << ELos1/ELos2 -1
	   <<  "  go(pbar)/bragg(p)-1= " << ELos3/ELos1 -1 << "   " 
	   << mnames[j] << "  Bethe(MeV*cm^2/g)= " 
	   << ELos2*g/(cm2*MeV*mat->GetDensity()) << G4endl;    
  }
  Energy = 10.*MeV;
  G4cout << "### E(MeV)= " << Energy/MeV << "  for " << pName 
	 << std::setprecision(5) << G4endl;
  for(G4int j=0; j<98; ++j) {
    const G4Material* mat = nist->FindOrBuildMaterial(mnames[j]);
    G4double Z = mat->GetZ();
    ELos1   = bragg->ComputeDEDXPerVolume(mat,prot,Energy,Ecut);     
    ELos2   = bethe1->ComputeDEDXPerVolume(mat,prot,Energy,Ecut);     
    ELos3   = qo->ComputeDEDXPerVolume(mat,part,Energy,Ecut);
    G4cout << "Z= " << Z << " p bragg/bethe-1= " << ELos1/ELos2 -1
	   <<  "  go(pbar)/bragg(p)-1= " << ELos3/ELos1 -1 << "   " 
	   << mnames[j] << "  Bethe(MeV*cm^2/g)= " 
	   << ELos2*g/(cm2*MeV*mat->GetDensity()) << G4endl;    
  }

  part = G4AntiProton::AntiProton();
  pName = part->GetParticleName();
  Energy = 2.*MeV;
  G4cout << "### E(MeV)= " << Energy/MeV << "  for " << pName 
	 << std::setprecision(5) << G4endl;
  for(G4int j=0; j<98; ++j) {
    const G4Material* mat = nist->FindOrBuildMaterial(mnames[j]);
    G4double Z = mat->GetZ();
    ELos1   = bragg->ComputeDEDXPerVolume(mat,prot,Energy,Ecut);     
    ELos2   = bethe2->ComputeDEDXPerVolume(mat,part,Energy,Ecut);     
    ELos3   = qo->ComputeDEDXPerVolume(mat,part,Energy,Ecut);
    G4cout << "Z= " << Z << "  bragg(p)/bethe(pbar)-1= " << ELos1/ELos2 -1
	   <<  " pbar go/bethe-1= " << ELos3/ELos2 -1 << "   " 
	   << mnames[j] << "  Bethe(MeV*cm^2/g)= " 
	   << ELos2*g/(cm2*MeV*mat->GetDensity()) << G4endl;    
  }
  Energy = 10.*MeV;
  G4cout << "### E(MeV)= " << Energy/MeV << "  for " << pName  
	 << std::setprecision(5) << G4endl;
  for(G4int j=0; j<98; ++j) {
    const G4Material* mat = nist->FindOrBuildMaterial(mnames[j]);
    G4double Z = mat->GetZ();
    ELos1   = bragg->ComputeDEDXPerVolume(mat,prot,Energy,Ecut);     
    ELos2   = bethe2->ComputeDEDXPerVolume(mat,part,Energy,Ecut);     
    ELos3   = qo->ComputeDEDXPerVolume(mat,part,Energy,Ecut);
    G4cout << "Z= " << Z << "  bragg(p)/bethe(pbar)-1= " << ELos1/ELos2 -1
	   <<  " pbar go/bethe-1= " << ELos3/ELos2 -1 << "   " 
	   << mnames[j] << "  Bethe(MeV*cm^2/g)= " 
	   << ELos2*g/(cm2*MeV*mat->GetDensity()) << G4endl;    
  }

  // Write to file
  for (size_t i = 0; i < nBinTab+1; i++) {
    Energy = pVector->GetLowEdgeEnergy(i); 
    ELos1   = bragg->ComputeDEDXPerVolume(material,part,Energy,Ecut);     
    ELos2   = bethe2->ComputeDEDXPerVolume(material,part,Energy,Ecut);     
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

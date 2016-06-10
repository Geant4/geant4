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
// $Id: G4EmSaturation.cc 81936 2014-06-06 15:42:55Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4EmSaturation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 18.02.2008
//
// Modifications:
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmSaturation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmSaturation::G4EmSaturation(G4int verb) 
  : manager(0)
{
  verbose = verb;
  manager = 0;

  curMaterial = 0;
  curBirks    = 0.0;
  curRatio    = 1.0;
  curChargeSq = 1.0;
  nMaterials  = nWarnings = 0;

  electron = 0;
  proton   = 0;
  nist     = G4NistManager::Instance();

  InitialiseG4materials(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmSaturation::~G4EmSaturation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmSaturation::VisibleEnergyDeposition(
                                      const G4ParticleDefinition* p, 
				      const G4MaterialCutsCouple* couple, 
				      G4double length,
				      G4double edep,
				      G4double niel)
{
  if(edep <= 0.0) { return 0.0; }

  G4double evis = edep;
  G4double bfactor = FindBirksCoefficient(couple->GetMaterial());

  if(bfactor > 0.0) { 

    G4int pdgCode = p->GetPDGEncoding();
    // atomic relaxations for gamma incident
    if(22 == pdgCode && electron) {
      //G4cout << "%% gamma edep= " << edep/keV << " keV " <<manager << G4endl; 
      evis /= (1.0 + bfactor*edep/manager->GetRange(electron,edep,couple));

      // energy loss
    } else {

      // protections
      G4double nloss = niel;
      if(nloss < 0.0) { nloss = 0.0; }
      G4double eloss = edep - nloss;

      // neutrons and neutral hadrons
      if(0.0 == p->GetPDGCharge() || eloss < 0.0 || length <= 0.0) {
	nloss = edep;
        eloss = 0.0;
      }

      // continues energy loss
      if(eloss > 0.0) { eloss /= (1.0 + bfactor*eloss/length); }
 
      // non-ionizing energy loss
      if(nloss > 0.0 && proton) {
        G4double escaled = nloss*curRatio;
	/*
        G4cout << "%% p edep= " << nloss/keV << " keV  Escaled= " 
	       << escaled << " MeV  in " << couple->GetMaterial()->GetName()
	       << "  " << p->GetParticleName()
	       << G4endl; 
	*/
        G4double range = manager->GetRange(proton,escaled,couple)/curChargeSq; 
	nloss /= (1.0 + bfactor*nloss/range);
      }

      evis = eloss + nloss;
    }
  }
  
  return evis;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmSaturation::FindG4BirksCoefficient(const G4Material* mat)
{
  G4String name = mat->GetName();
  // is this material in the vector?
  
  for(G4int j=0; j<nG4Birks; ++j) {
    if(name == g4MatNames[j]) {
      if(verbose > 0) 
	G4cout << "### G4EmSaturation::FindG4BirksCoefficient for "
	       << name << " is " << g4MatData[j]*MeV/mm << " mm/MeV "
	       << G4endl;
      return g4MatData[j];
    }
  }
  return FindBirksCoefficient(mat);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmSaturation::InitialiseBirksCoefficient(const G4Material* mat)
{
  // electron and proton should exist in any case
  if(!manager) {
    manager = G4LossTableManager::Instance();
    electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  }

  curMaterial = mat;
  curBirks = 0.0;
  curRatio = 1.0;
  curChargeSq = 1.0;

  // seach in the run-time list
  for(G4int i=0; i<nMaterials; ++i) {
    if(mat == matPointers[i]) {
      curBirks = mat->GetIonisation()->GetBirksConstant();
      curRatio = massFactors[i];
      curChargeSq = effCharges[i];
      return;
    }
  }

  G4String name = mat->GetName();
  curBirks = mat->GetIonisation()->GetBirksConstant();

  // material has no Birks coeffitient defined
  // seach in the Geant4 list
  if(curBirks == 0.0) {
    for(G4int j=0; j<nG4Birks; ++j) {
      if(name == g4MatNames[j]) {
	mat->GetIonisation()->SetBirksConstant(g4MatData[j]);
        curBirks = g4MatData[j];
        break;
      }
    }
  }

  if(curBirks == 0.0) {
    if(0 < nWarnings) {
      ++nWarnings;
      G4ExceptionDescription ed;
      ed << "Birks constants are not defined for material " << name 
	 << " ! \n Define Birks constants for the material" 
	 << " or not apply saturation.";
      G4Exception("G4EmSaturation::InitialiseBirksCoefficient", "em0088",
		  JustWarning, ed);
    }
    return;
  }

  // compute mean mass ratio
  curRatio = 0.0;
  curChargeSq = 0.0;
  G4double norm = 0.0;
  const G4ElementVector* theElementVector = mat->GetElementVector();
  const G4double* theAtomNumDensityVector = mat->GetVecNbOfAtomsPerVolume();
  size_t nelm = mat->GetNumberOfElements();
  for (size_t i=0; i<nelm; ++i) {
    const G4Element* elm = (*theElementVector)[i];
    G4double Z = elm->GetZ();
    G4double w = Z*Z*theAtomNumDensityVector[i];
    curRatio += w/nist->GetAtomicMassAmu(G4int(Z));
    curChargeSq = Z*Z*w;
    norm += w;
  }
  curRatio *= proton_mass_c2/norm;
  curChargeSq /= norm;

  // store results
  matPointers.push_back(mat);
  matNames.push_back(name);
  massFactors.push_back(curRatio);
  effCharges.push_back(curChargeSq);
  nMaterials++;
  if(verbose > 0) {
    G4cout << "### G4EmSaturation::FindBirksCoefficient Birks coefficient for "
	   << name << "  " << curBirks*MeV/mm << " mm/MeV" << G4endl;
  }
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmSaturation::DumpBirksCoefficients()
{
  if(nMaterials > 0) {
    G4cout << "### Birks coeffitients used in run time" << G4endl;
    for(G4int i=0; i<nMaterials; ++i) {
      G4double br = matPointers[i]->GetIonisation()->GetBirksConstant();
      G4cout << "   " << matNames[i] << "     " 
	     << br*MeV/mm << " mm/MeV" << "     "
	     << br*matPointers[i]->GetDensity()*MeV*cm2/g 
	     << " g/cm^2/MeV" 
	     << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmSaturation::DumpG4BirksCoefficients()
{
  if(nG4Birks > 0) {
    G4cout << "### Birks coeffitients for Geant4 materials" << G4endl;
    for(G4int i=0; i<nG4Birks; ++i) {
      G4cout << "   " << g4MatNames[i] << "   " 
	     << g4MatData[i]*MeV/mm << " mm/MeV" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmSaturation::InitialiseG4materials()
{
  // M.Hirschberg et al., IEEE Trans. Nuc. Sci. 39 (1992) 511
  // SCSN-38 kB = 0.00842 g/cm^2/MeV; rho = 1.06 g/cm^3
  g4MatNames.push_back("G4_POLYSTYRENE");
  g4MatData.push_back(0.07943*mm/MeV);

  // C.Fabjan (private communication)
  // kB = 0.006 g/cm^2/MeV; rho = 7.13 g/cm^3
  g4MatNames.push_back("G4_BGO");
  g4MatData.push_back(0.008415*mm/MeV);

  // A.Ribon analysis of publications
  // Scallettar et al., Phys. Rev. A25 (1982) 2419.
  // NIM A 523 (2004) 275. 
  // kB = 0.022 g/cm^2/MeV; rho = 1.396 g/cm^3; 
  // ATLAS Efield = 10 kV/cm provide the strongest effect 
  g4MatNames.push_back("G4_lAr");
  g4MatData.push_back(0.1576*mm/MeV);

  //G4_BARIUM_FLUORIDE
  //G4_CESIUM_IODIDE
  //G4_GEL_PHOTO_EMULSION
  //G4_PHOTO_EMULSION
  //G4_PLASTIC_SC_VINYLTOLUENE
  //G4_SODIUM_IODIDE
  //G4_STILBENE
  //G4_lAr
  //G4_PbWO4
  //G4_Lucite

  nG4Birks = g4MatData.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

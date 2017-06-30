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
// $Id: G4EmSaturation.cc 104372 2017-05-29 09:55:44Z gcosmo $
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

G4int G4EmSaturation::nMaterials = 0;
std::vector<G4double> G4EmSaturation::massFactors;
std::vector<G4double> G4EmSaturation::effCharges;
std::vector<G4double> G4EmSaturation::g4MatData;
std::vector<G4String> G4EmSaturation::g4MatNames;

G4EmSaturation::G4EmSaturation(G4int verb) 
{
  verbose = verb;

  nWarnings = nG4Birks = 0;

  electron = nullptr;
  proton   = nullptr;
  nist     = G4NistManager::Instance();
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
                                      G4double niel) const
{
  if(edep <= 0.0) { return 0.0; }

  G4double evis = edep;
  G4double bfactor = couple->GetMaterial()->GetIonisation()->GetBirksConstant();

  if(bfactor > 0.0) { 

    // atomic relaxations for gamma incident
    if(22 ==  p->GetPDGEncoding()) {
      //G4cout << "%% gamma edep= " << edep/keV << " keV " <<manager << G4endl; 
      evis /= (1.0 + bfactor*edep/
        G4LossTableManager::Instance()->GetRange(electron,edep,couple));

      // energy loss
    } else {

      // protections
      G4double nloss = std::max(niel, 0.0);
      G4double eloss = edep - nloss;

      // neutrons and neutral hadrons
      if(0.0 == p->GetPDGCharge() || eloss < 0.0 || length <= 0.0) {
        nloss = edep;
        eloss = 0.0;
      } else {

	// continues energy loss
	eloss /= (1.0 + bfactor*eloss/length); 
      }
      // non-ionizing energy loss
      if(nloss > 0.0) {
        G4int idx = couple->GetMaterial()->GetIndex();
        G4double escaled = nloss*massFactors[idx];
        /*
        G4cout << "%% p edep= " << nloss/keV << " keV  Escaled= " 
               << escaled << " MeV  in " << couple->GetMaterial()->GetName()
               << "  " << p->GetParticleName()
               << G4endl; 
        */
        G4double range = G4LossTableManager::Instance()
          ->GetRange(proton,escaled,couple)/effCharges[idx]; 
        nloss /= (1.0 + bfactor*nloss/range);
      }
      evis = eloss + nloss;
    }
  }
  return evis;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmSaturation::InitialiseG4Saturation()
{
  nMaterials = G4Material::GetNumberOfMaterials();
  massFactors.resize(nMaterials, 1.0);
  effCharges.resize(nMaterials, 1.0);

  if(0 == nG4Birks) {  InitialiseG4materials(); }

  for(G4int i=0; i<nMaterials; ++i) {
    InitialiseBirksCoefficient((*G4Material::GetMaterialTable())[i]);
  }
  if(verbose > 0) { DumpBirksCoefficients(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmSaturation::FindG4BirksCoefficient(const G4Material* mat)
{
  if(0 == nG4Birks) {  InitialiseG4materials(); }

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
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmSaturation::InitialiseBirksCoefficient(const G4Material* mat)
{
  // electron and proton should exist in any case
  if(!electron) {
    electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
    if(!electron || !proton) {
      G4Exception("G4EmSaturation::InitialiseBirksCoefficient", "em0001",
		  FatalException, "both electron and proton should exist");
    }
  }

  G4double curBirks = mat->GetIonisation()->GetBirksConstant();

  G4String name = mat->GetName();

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

  if(curBirks == 0.0) { return; }

  // compute mean mass ratio
  G4double curRatio = 0.0;
  G4double curChargeSq = 0.0;
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
  G4int idx = mat->GetIndex();
  massFactors[idx] = curRatio;
  effCharges[idx] = curChargeSq;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmSaturation::DumpBirksCoefficients()
{
  G4cout << "### Birks coefficients used in run time" << G4endl;
  const G4MaterialTable* mtable = G4Material::GetMaterialTable();
  for(G4int i=0; i<nMaterials; ++i) {
    const G4Material* mat = (*mtable)[i];
    G4double br = mat->GetIonisation()->GetBirksConstant();
    if(br > 0.0) {
      G4cout << "   " << mat->GetName() << "     " 
	     << br*MeV/mm << " mm/MeV" << "     "
	     << br*mat->GetDensity()*MeV*cm2/g 
	     << " g/cm^2/MeV  massFactor=  " << massFactors[i]  
	     << " effCharge= " << effCharges[i] << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmSaturation::DumpG4BirksCoefficients()
{
  if(nG4Birks > 0) {
    G4cout << "### Birks coefficients for Geant4 materials" << G4endl;
    for(G4int i=0; i<nG4Birks; ++i) {
      G4cout << "   " << g4MatNames[i] << "   " 
             << g4MatData[i]*MeV/mm << " mm/MeV" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmSaturation::InitialiseG4materials()
{
  nG4Birks = 4;
  g4MatData.reserve(nG4Birks);

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
  // kB = 0.1576*mm/MeV
  // A. Kiryunin and P.Strizenec "Geant4 hadronic 
  // working group meeting " kB = 0.041/9.13 g/cm^2/MeV  
  g4MatNames.push_back("G4_lAr");
  g4MatData.push_back(0.032*mm/MeV);

  //G4_BARIUM_FLUORIDE
  //G4_CESIUM_IODIDE
  //G4_GEL_PHOTO_EMULSION
  //G4_PHOTO_EMULSION
  //G4_PLASTIC_SC_VINYLTOLUENE
  //G4_SODIUM_IODIDE
  //G4_STILBENE
  //G4_lAr

  //G4_PbWO4 - CMS value
  g4MatNames.push_back("G4_PbWO4");
  g4MatData.push_back(0.0333333*mm/MeV);

  //G4_Lucite

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

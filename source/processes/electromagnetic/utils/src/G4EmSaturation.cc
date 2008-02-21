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
// $Id: G4EmSaturation.cc,v 1.4 2008-02-21 18:11:19 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4LossTableManager.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmSaturation::G4EmSaturation()
{
  verbose = 1;
  pointers = false;
  curMaterial = 0;
  curBirks = 0.0;
  nMaterials = 0;
  Initialise();
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
  if(edep <= 0.0) return 0.0;

  G4double evis = edep;
  G4double bfactor = FindBirksCoefficient(couple->GetMaterial());

  if(bfactor > 0.0) { 

    if(!pointers) {
      pointers = true;
      manager = G4LossTableManager::Instance();
      gamma = G4Gamma::Gamma();
      electron = G4Electron::Electron();
      proton = G4Proton::Proton();
      neutron = G4Neutron::Neutron();
    }

    // atomic relaxations
    if(p == gamma) {
      evis /= (1.0 + bfactor*edep/manager->GetRange(electron,edep,couple));

      // energy loss
    } else {

      // protections
      G4double nloss = niel;
      if(nloss < 0.0) nloss = 0.0;
      G4double eloss = edep - nloss;
      if(p == neutron || eloss < 0.0 || length <= 0.0) {
	nloss = edep;
        eloss = 0.0;
      }

      // continues energy loss
      if(eloss > 0.0) eloss /= (1.0 + bfactor*eloss/length);
 
      // non-ionizing energy loss
      if(nloss > 0.0) {
        G4double q = p->GetPDGCharge()/eplus;
        G4double escaled = nloss*proton_mass_c2/curAtomicMass;
        G4double s = manager->GetRange(proton,escaled,couple)/(q*q); 
	nloss /= (1.0 + bfactor*nloss/s);
      }

      evis = eloss + nloss;
    }
  }
  
  return evis;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmSaturation::FillBirksCoefficient(const G4Material* mat, G4double birks)
{
  G4String name = mat->GetName();
  // is this material in the vector?
  for(G4int i=0; i<nMaterials; i++) {
    if(mat == matPointers[i]) {
      if(verbose > 0) 
	G4cout << "### G4EmSaturation::FillBirksCoefficient Birks coefficient for "
	       << name << " BirksOld= " << matData[i] 
	       << " is substituted by " << birks*MeV/mm << " mm/MeV" << G4endl;
      matData[i] = birks;
      return;
    }
  }
  matPointers.push_back(mat);
  matNames.push_back(name);
  matData.push_back(birks);
  nMaterials++;
  if(verbose > 0) 
    G4cout << "### G4EmSaturation::FillBirksCoefficient Birks coefficient for "
	   << name << "  " << birks*MeV/mm << " mm/MeV" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmSaturation::FindG4BirksCoefficient(const G4Material* mat)
{
  G4String name = mat->GetName();
  // is this material in the vector?
  
  for(G4int i=0; i<nMaterials; i++) {
    if(mat == matPointers[i] || name == matNames[i]) {
      if(verbose > 0) 
	G4cout << "### G4EmSaturation::FindG4BirksCoefficient for "
	       << name << " find out " << matPointers[i]->GetName() 
	       << " with coefficient " << matData[i] << G4endl;
      return matData[i];
    }
  }
  return FindBirksCoefficient(mat);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmSaturation::FindBirksCoefficient(const G4Material* mat)
{
  if(mat == curMaterial) return curBirks;

  curMaterial = mat;
  curBirks = 0.0;

  // seach in the run-time list
  for(G4int i=0; i<nMaterials; i++) {
    if(mat == matPointers[i]) {
      curBirks = matData[i];
      return curBirks;
    }
  }

  // seach in the Geant4 list
  G4String name = mat->GetName();
  for(G4int j=0; j<nG4Birks; j++) {
    if(name == g4MatNames[j]) {
      curBirks = g4MatData[j];
      matPointers.push_back(mat);
      matNames.push_back(name);
      matData.push_back(curBirks);
      nMaterials++;
      if(verbose > 0) 
	G4cout << "### G4EmSaturation::FindBirksCoefficient Birks coefficient for "
	       << name << "  " << curBirks*MeV/mm << " mm/MeV" << G4endl;
      return curBirks;
    }
  }

  if(verbose > 0) 
    G4cout << "### G4EmSaturation::FindBirksCoefficient Birks coefficient fails for "
	   << name << G4endl;
  return curBirks;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmSaturation::DumpBirksCoefficients()
{
  if(nMaterials > 0) {
    G4cout << "### Birks coeffitients used in run time" << G4endl;
    for(G4int i=0; i<nMaterials; i++) {
      G4cout << "   " << matNames[i] << "     " 
	     << matData[i]*MeV/mm << " mm/MeV" << "     "
	     << matData[i]*matPointers[i]->GetDensity()*MeV*cm2/g << " g/cm^2/MeV" 
	     << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmSaturation::DumpG4BirksCoefficients()
{
  if(nG4Birks > 0) {
    G4cout << "### Birks coeffitients for Geant4 materials" << G4endl;
    for(G4int i=0; i<nG4Birks; i++) {
      G4cout << "   " << g4MatNames[i] << "   " 
	     << g4MatData[i]*MeV/mm << " mm/MeV" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmSaturation::Initialise()
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

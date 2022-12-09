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
//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Description: Data on stopping power
//
// Author:      Vladimir Ivanchenko
//
// Creation date: 23.10.2021
// 
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <fstream>
#include <sstream>

#include "G4IonICRU73Data.hh" 
#include "G4SystemOfUnits.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4EmParameters.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"
#include "G4Element.hh"

namespace {

const G4String namesICRU73[31] = {
"G4_A-150_TISSUE"
"G4_ADIPOSE_TISSUE_ICRP"
"G4_AIR"
"G4_ALUMINUM_OXIDE"
"G4_BONE_COMPACT_ICRU"
"G4_BONE_CORTICAL_ICRP"
"G4_C-552"
"G4_CALCIUM_FLUORIDE"
"G4_CARBON_DIOXIDE"
"G4_KAPTON"
"G4_LITHIUM_FLUORIDE"
"G4_LITHIUM_TETRABORATE"
"G4_LUCITE"
"G4_METHANE"
"G4_MUSCLE_STRIATED_ICRU"
"G4_MYLAR"
"G4_NYLON-6-6"
"G4_PHOTO_EMULSION"
"G4_PLASTIC_SC_VINYLTOLUENE"
"G4_POLYCARBONATE"
"G4_POLYETHYLENE"
"G4_POLYSTYRENE"
"G4_PROPANE"
"G4_Pyrex_Glass"
"G4_SILICON_DIOXIDE"
"G4_SODIUM_IODIDE"
"G4_TEFLON"
"G4_TISSUE-METHANE"
"G4_TISSUE-PROPANE"
"G4_WATER"
"G4_WATER_VAPOR"};
  const G4String namesICRU90[3] = {"G4_AIR","G4_GRAPHITE","G4_WATER"};
  const G4double densityCoef[3] = {0.996, 1.025, 0.998};
  const G4int NZ = 27;
  const G4int zdat[NZ] = { 5,   6,  7,  8, 13, 14, 15, 16, 22, 26,
                          28, 29, 30, 31, 32, 33, 34, 40, 47, 48,
                          49, 51, 52, 72, 73, 74, 79 }; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4IonICRU73Data::G4IonICRU73Data()
{
  fEmin = 0.025*CLHEP::MeV;
  fEmax = 2.5*CLHEP::MeV;
  fNbins = fNbinsPerDecade*G4lrint(std::log10(fEmax/fEmin));
  fVector = new G4PhysicsFreeVector(fSpline);
  for(G4int i=0; i<81; ++i) {
    fMatData[i] = new std::vector<G4PhysicsLogVector*>;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4IonICRU73Data::~G4IonICRU73Data()
{
  delete fVector;
  for(G4int i=0; i<81; ++i) {
    auto v = fMatData[i];
    for(G4int j=0; j<fNmat; ++j) {
      delete (*v)[j];
    }
    delete v;
    for(G4int j=0; j<93; ++j) { delete fElmData[i][j]; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4IonICRU73Data::GetDEDX(const G4Material* mat, const G4int Z,
                                  const G4double e, const G4double loge) const
{
  G4PhysicsLogVector* v = nullptr;
  G4int Z2 = std::min(Z, 80);
  if(1 == mat->GetNumberOfElements()) {
    G4int Z1 = std::min((*(mat->GetElementVector()))[0]->GetZasInt(), 80);
    v = fElmData[Z2][Z1];
  } else {
    G4int idx = fMatIndex[mat->GetIndex()];
    if(idx < fNmat) {
      v = (*(fMatData[Z2]))[idx];
    }
  }
  if(nullptr == v) { return 0.0; }
  G4double res = (e > fEmin) ? v->LogVectorValue(e, loge) 
    : (*v)[0]*std::sqrt(e/fEmin);
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void G4IonICRU73Data::Initialise()
{
  // fill directory path
  if(fDataDirectory.empty()) {
    const char* path = G4FindDataDir("G4LEDATA");
    if (nullptr != path) {
      std::ostringstream ost;
      ost << path << "/ion_stopping_data/";
      fDataDirectory = ost.str();
    } else {
      G4Exception("G4IonICRU73Data::Initialise(..)","em013",
                  FatalException,
                  "Environment variable G4LEDATA is not defined");
    }
  }

  std::size_t nmat = G4Material::GetNumberOfMaterials();
  if(nmat == fMatIndex.size()) { return; }
 
  if(0 < fVerbose) {
    G4cout << "### G4IonICRU73Data::Initialise() for " << nmat
           << " materials" << G4endl;
  } 
  fMatIndex.resize(nmat, -1);
  for(G4int j=0; j<81; ++j) {
    fMatData[j]->resize(nmat, nullptr);
  }
  G4bool useICRU90 = G4EmParameters::Instance()->UseICRU90Data();
  auto mtable = G4Material::GetMaterialTable();

  for(G4int i=0; i<(G4int)nmat; ++i) {
    fNmat = i;
    const G4Material* mat = (*mtable)[i];
    G4int idx = (G4int)mat->GetIndex();
    if(1 < fVerbose) {
      G4cout << i << ".  material:" << mat->GetName() 
             << "  idx=" << idx << G4endl;
    }
    if(fMatIndex[idx] == -1) {
      fMatIndex[idx] = i;
      G4String matname = mat->GetName();
      G4bool isOK = false; 
      if(1 == mat->GetNumberOfElements()) {
        ReadElementData(mat, useICRU90);
        isOK = true;
        if(1 < fVerbose) {
          G4cout << "Material from single element" << G4endl;
        }
      }
      if(!isOK && useICRU90) {
        for(G4int j=0; j<3; ++j) {
          if(matname == namesICRU90[j]) {
            ReadMaterialData(mat, densityCoef[j], true);
            isOK = true;
            if(1 < fVerbose) {
              G4cout << "ICRU90 material" << G4endl;
            }
            break;
          }
        }
      }
      if(!isOK) {
        for(G4int j=0; j<31; ++j) {
          if(matname == namesICRU73[j]) {
            ReadMaterialData(mat, 1.0, false);
            isOK = true;
            if(1 < fVerbose) {
              G4cout << "ICRU73 material" << G4endl;
            }
            break;
          }
        }
      }
      if(!isOK) {
        ReadElementData(mat, useICRU90);
        if(1 < fVerbose) {
          G4cout << "Read element data" << G4endl;
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void G4IonICRU73Data::ReadMaterialData(const G4Material* mat, 
                                       const G4double coeff, 
                                       const G4bool useICRU90)
{
  G4String name = mat->GetName();
  for(G4int Z=3; Z<81; ++Z) {
    std::ostringstream ost;
    ost << fDataDirectory << "icru";
    G4int Z1 = Z;
    G4double scale = 1.0;
    if(useICRU90 && Z <= 18) {
      ost << "90";
    } else {
      ost << "73";
      for(G4int i=0; i<NZ; ++i) {
	if(Z == zdat[i]) {
	  break;
	} else if(i == NZ-1) {
	  Z1 = zdat[NZ - 1];
	  scale = (G4double)(Z*Z)/(G4double)(Z1*Z1);
	} else if(Z > zdat[i] && Z < zdat[i+1]) {
	  if(Z - zdat[i] <= zdat[i + 1] - Z) {
	    Z1 = zdat[i];
	  } else {
	    Z1 = zdat[i+1];
	  }
	  scale = (G4double)(Z*Z)/(G4double)(Z1*Z1);
	  break;
	}
      }
    }
    if(nullptr == (*(fMatData[Z1]))[fNmat]) {
      ost << "/z" << Z1 << "_" << name << ".dat"; 
      G4PhysicsLogVector* v = RetrieveVector(ost, false);
      if(nullptr != v) {
	const G4double fact = coeff *
	  mat->GetDensity() * CLHEP::MeV * 1000 * CLHEP::cm2 / CLHEP::g; 
	v->ScaleVector(CLHEP::MeV, fact);
	if(2 < fVerbose) {
	  G4cout << "### Data for " << name 
		 << " and projectile Z=" << Z1 << G4endl;
	  G4cout << *v << G4endl;
	}
	(*(fMatData[Z1]))[fNmat] = v;
      }
    }
    if(Z != Z1) {
      auto v2 = (*(fMatData[Z1]))[fNmat];
      if(nullptr != v2) {
	auto v1 = new G4PhysicsLogVector(*v2);
	(*(fMatData[Z]))[fNmat] = v1;
	v1->ScaleVector(1.0, scale);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void G4IonICRU73Data::ReadElementData(const G4Material* mat, G4bool useICRU90)
{
  const G4ElementVector* elmv = mat->GetElementVector();
  const G4double* dens = mat->GetFractionVector();
  const G4int nelm = (G4int)mat->GetNumberOfElements();
  for(G4int Z=3; Z<81; ++Z) {
    G4PhysicsLogVector* v = nullptr; 
    if(1 == nelm) {
      v = FindOrBuildElementData(Z, (*elmv)[0]->GetZasInt(), useICRU90);
    } else {
      v = new G4PhysicsLogVector(fEmin, fEmax, fNbins, fSpline);
      for(G4int i=0; i<=fNbins; ++i) {
        G4double dedx = 0.0;
        for(G4int j=0; j<nelm; ++j) {
          G4PhysicsLogVector* v1 = 
            FindOrBuildElementData(Z, (*elmv)[j]->GetZasInt(), useICRU90);
          dedx += (*v1)[i]*dens[j];
        }
        v->PutValue(i, dedx);
      }
      if(fSpline) { v->FillSecondDerivatives(); }
    }
    (*(fMatData[Z]))[fNmat] = v;
    // scale data for correct units
    if(nullptr != v) {
      const G4double fact =
        mat->GetDensity() * CLHEP::MeV * 1000 * CLHEP::cm2 / CLHEP::g; 
      v->ScaleVector(CLHEP::MeV, fact);
      if(2 < fVerbose) {
        G4cout << "### Data for "<< mat->GetName() 
               << " for projectile Z=" << Z << G4endl;
        G4cout << *v << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4PhysicsLogVector* 
G4IonICRU73Data::FindOrBuildElementData(const G4int Z, const G4int Z1, 
                                        G4bool useICRU90)
{
  G4PhysicsLogVector* v = nullptr;
  if(Z <= 80 && Z1 <= 92) {
    v = fElmData[Z][Z1];
    if(nullptr == v) {
      G4int Z2 = Z1;
      G4bool isICRU90 = (useICRU90 && Z <= 18 && 
                         (Z1 == 1 || Z1 == 6 || Z1 == 7 || Z1 == 8));

      G4double scale = 1.0;      
      if(!isICRU90) {
        for(G4int i=0; i<NZ; ++i) {
          if(Z1 == zdat[i]) {
	    break;
	  } else if(i == NZ-1) {
	    Z2 = zdat[NZ - 1];
            scale = (G4double)Z1/(G4double)Z2;
	  } else if(Z1 > zdat[i] && Z1 < zdat[i+1]) {
            if(Z1 - zdat[i] <= zdat[i + 1] - Z1) {
	      Z2 = zdat[i];
	    } else {
	      Z2 = zdat[i+1];
	    }
            scale = (G4double)Z1/(G4double)Z2;
	    break;
	  }
	}
      }

      std::ostringstream ost;
      ost << fDataDirectory << "icru";
      if(isICRU90) { ost << "90"; }
      else { ost << "73"; }
      ost << "/z" << Z << "_" << Z2 << ".dat"; 
      v = RetrieveVector(ost, false);
      fElmData[Z][Z2] = v;
      if(Z1 != Z2 && nullptr != v) {
	fElmData[Z][Z1] = new G4PhysicsLogVector(*v);
	fElmData[Z][Z1]->ScaleVector(1.0, scale);
      }
    }
  }
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4PhysicsLogVector*
G4IonICRU73Data::RetrieveVector(std::ostringstream& ost, G4bool warn)
{
  G4PhysicsLogVector* v = nullptr;
  std::ifstream filein(ost.str().c_str());
  if (!filein.is_open()) {
    if(warn) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
         << "> is not opened";
      G4Exception("G4IonICRU73Data::RetrieveVector(..)","em013",
                  FatalException, ed, "Check G4LEDATA");
    }
  } else {
    if(fVerbose > 0) {
      G4cout << "File " << ost.str() 
             << " is opened by G4IonICRU73Data" << G4endl;
    }

    // retrieve data from DB
    if(!fVector->Retrieve(filein, true)) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
         << "> is not retrieved!";
      G4Exception("G4IonICRU73Data::RetrieveVector(..)","had015",
                  FatalException, ed, "Check G4LEDATA");
    } else {
      if(fSpline) { fVector->FillSecondDerivatives(); }
      v = new G4PhysicsLogVector(fEmin, fEmax, fNbins, fSpline);
      for(G4int i=0; i<=fNbins; ++i) {
        G4double e = v->Energy(i);
        G4double dedx = fVector->Value(e);
        v->PutValue(i, dedx);
      }
      if(fSpline) { v->FillSecondDerivatives(); }
      if(fVerbose > 1) { G4cout << *v << G4endl; }
    }
  }
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

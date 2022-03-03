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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4IonICRU73Data::G4IonICRU73Data()
{
  fEmin = 0.025*CLHEP::MeV;
  fEmax = 2.5*CLHEP::MeV;
  fNbins = fNbinsPerDecade*G4lrint(std::log10(fEmax/fEmin) + 0.000001);
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
    for(G4int j=0; j<81; ++j) { delete fElmData[i][j]; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4IonICRU73Data::GetDEDX(const G4Material* mat, const G4int Z,
				  const G4double e, const G4double loge) const
{
  G4PhysicsLogVector* v = nullptr;
  if(1 == mat->GetNumberOfElements()) {
    G4int Z1 = (*(mat->GetElementVector()))[0]->GetZasInt();
    if(Z <= 80 && Z1 <= 80) {
      v = fElmData[Z][Z1];
    }
  } else {
    G4int idx = fMatIndex[mat->GetIndex()];
    if(idx < fNmat && Z <= 80) {
      v = (*(fMatData[Z]))[idx];
    }
  }
  G4double res = (nullptr != v) ? v->LogVectorValue(e, loge) : 0.0;
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void G4IonICRU73Data::Initialise()
{
  // fill directory path
  if(fDataDirectory.empty()) {
    char* path = std::getenv("G4LEDATA");
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

  G4bool useICRU90 = G4EmParameters::Instance()->UseICRU90Data();
  size_t nmat = G4Material::GetNumberOfMaterials();
  fMatIndex.resize(nmat, -1);

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  for(size_t i=0; i<numOfCouples; ++i) {
    const G4MaterialCutsCouple* couple = 
      theCoupleTable->GetMaterialCutsCouple(i);
    const G4Material* mat = couple->GetMaterial();
    G4int idx = mat->GetIndex();
    if(fMatIndex[idx] == -1) {
      G4String matname = mat->GetName();
      G4bool isOK = false; 
      if(1 == mat->GetNumberOfElements()) {
	ReadElementData(mat, useICRU90);
        isOK = true;
      }
      if(!isOK && useICRU90) {
	for(G4int j=0; j<3; ++j) {
	  if(matname == namesICRU90[j]) {
	    ReadMaterialData(matname, true);
            isOK = true;
	    fMatIndex[idx] = fNmat; 
	    ++fNmat;
	    break;
	  }
	}
      }
      if(!isOK) {
	for(G4int j=0; j<31; ++j) {
	  if(matname == namesICRU73[j]) {
	    ReadMaterialData(matname, false);
            isOK = true;
	    fMatIndex[idx] = fNmat; 
	    ++fNmat;
	    break;
	  }
	}
      }
      if(!isOK) {
	ReadElementData(mat, useICRU90);
	fMatIndex[idx] = fNmat; 
	++fNmat;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void G4IonICRU73Data::ReadMaterialData(const G4String& name, G4bool isICRU90)
{
  for(G4int Z=3; Z<81; ++Z) {
    std::ostringstream ost;
    ost << fDataDirectory << "icru";
    if(isICRU90 && Z <= 18) { ost << "90"; }
    else { ost << "73"; }
    ost << "/z" << Z << "_" << name << ".dat"; 
    G4PhysicsLogVector* v = RetrieveVector(ost, false);
    (fMatData[Z])->push_back(v);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void G4IonICRU73Data::ReadElementData(const G4Material* mat, G4bool useICRU90)
{
  const G4ElementVector* elmv = mat->GetElementVector();
  const G4double* dens = mat->GetFractionVector();
  const G4int nelm = mat->GetNumberOfElements();
  for(G4int Z=3; Z<81; ++Z) {
    if(1 == nelm) {
      FindOrBuildElementData(Z, (*elmv)[0]->GetZasInt(), useICRU90);
    } else {
      G4PhysicsLogVector* v2 = nullptr; 
      for(G4int j=0; j<nelm; ++j) {
	v2 = FindOrBuildElementData(Z, (*elmv)[j]->GetZasInt(), useICRU90);
        if(nullptr == v2) { break; }
      }
      // for one or more elements there is no data
      if(nullptr == v2) {
	fMatData[Z]->push_back(v2);
        continue;
      }
      G4PhysicsLogVector* v =
	new G4PhysicsLogVector(fEmin, fEmax, fNbins, fSpline);
      for(G4int i=0; i<=fNbins; ++i) {
	G4double dedx = 0;
        for(G4int j=0; j<nelm; ++j) {
	  G4PhysicsLogVector* v1 = 
	    FindOrBuildElementData(Z, (*elmv)[j]->GetZasInt(), useICRU90);
          dedx += (*v1)[i]*dens[j];
	}
        v->PutValue(i, dedx);
      }
      if(fSpline) { v->FillSecondDerivatives(); }
      fMatData[Z]->push_back(v);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4PhysicsLogVector* 
G4IonICRU73Data::FindOrBuildElementData(const G4int Z, const G4int Z1, 
                                        G4bool useICRU90)
{
  G4PhysicsLogVector* v = nullptr;
  if(Z <= 80 && Z1 <= 80) {
    v = fElmData[Z][Z1];
    if(nullptr == v) {
      G4bool isICRU90 = (useICRU90 && Z <= 18 && 
			 (Z1 == 1 || Z1 == 6 || Z1 == 7 || Z1 == 8)); 
      std::ostringstream ost;
      ost << fDataDirectory << "icru";
      if(isICRU90) { ost << "90"; }
      else { ost << "73"; }
      ost << "/z" << Z << "_" << Z1 << ".dat"; 
      v = RetrieveVector(ost, false);
      fElmData[Z][Z1] = v;
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
      const G4double fact = CLHEP::MeV * CLHEP::cm2 /( 0.001 * CLHEP::g); 
      v->ScaleVector(CLHEP::MeV, fact);
      if(fSpline) { v->FillSecondDerivatives(); }
      if(fVerbose > 1) { G4cout << *v << G4endl; }
    }
  }
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
// $Id: G4eIonisationParameters.cc,v 1.3 2001-09-14 08:17:12 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
// 12.09.01 V.Ivanchenko    Add param and interpolation of parameters  
//
// -------------------------------------------------------------------

#include "G4eIonisationParameters.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4ShellEMDataSet.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4MaterialTable.hh"
#include "G4DataVector.hh"
#include "g4std/fstream"
#include "g4std/strstream"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eIonisationParameters:: G4eIonisationParameters(G4int minZ, G4int maxZ)
  : zMin(minZ), zMax(maxZ),
  length(16),
  interpolation(0)
{
  LoadData();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eIonisationParameters::~G4eIonisationParameters()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eIonisationParameters::Parameter(G4int Z, G4int shellIndex, 
                                                     G4int parameterIndex, 
                                                     G4double e) const
{
  G4double value = 0.;
  G4int id = Z*20 + parameterIndex;
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;

  pos = param.find(id);
  if (pos!= param.end()) {
    G4VEMDataSet* dataSet = pos->second;
    size_t nShells = dataSet->NumberOfComponents();

    if(shellIndex < nShells) { 
      const G4VEMDataSet* comp = dataSet->GetComponent(shellIndex);
      const G4DataVector ener = comp->GetEnergies(0);
      G4double ee = G4std::max(ener.front(),G4std::min(ener.back(),e));
      G4double value = comp->FindValue(ee);

    } else {
      G4cout << "WARNING: G4IonisationParameters::FindParameters "
             << "has no parameters for shell= " << shellIndex 
             << "; Z= " << Z
             << G4endl;
    }
  } else {
    G4cout << "WARNING: G4IonisationParameters::FindValue "
           << "did not find ID = "
           << index << G4endl;
  }

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4DataVector& G4eIonisationParameters::Parameters(G4int Z, 
                                                       G4int shellIndex,
                                                       G4double e) const
{
// To be implemented
  G4DataVector* p = new G4DataVector();
  return *p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eIonisationParameters::LoadData()
{
  if(!interpolation)  interpolation  = new G4SemiLogInterpolation();
  // define active elements

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
     G4Exception("G4CrossSectionHandler: no MaterialTable found)");

  G4int nMaterials = materialTable->length();
  
  for (G4int m=0; m<nMaterials; m++) {

    const G4Material* material= (*materialTable)[m];        
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4int nElements = material->GetNumberOfElements();
      
    for (size_t iEl=0; iEl<nElements; iEl++) {
      G4Element* element = (*elementVector)[iEl];
      G4double Z = element->GetZ();
      if (!(activeZ.contains(Z)) && Z > 0 && Z < zMax) {
	 activeZ.push_back(Z);
      }
    }
  }
  char* path = getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep = "G4EMDataSet - G4LEDATA environment variable not set";
      G4Exception(excep);
    }
  G4String pathString(path);
  pathString += "/ioni/io-co-";  

  G4double energy;
  G4DataVector* e;
  G4VEMDataSet* set;
  G4std::vector<G4DataVector*> a;
  G4std::vector<G4ShellEMDataSet*> p;
  a.resize(length);
  p.resize(length);

  size_t nZ = activeZ.size();

  for (G4int i=0; i<nZ; i++) {

    G4int Z = (G4int)activeZ[i];
    char nameChar[100] = {""};
    G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
    ost << pathString << Z << ".dat";
    G4String name(nameChar);  

    G4std::ifstream file(name);
    G4std::filebuf* lsdp = file.rdbuf();
  
    if (! (lsdp->is_open()) ) {
      G4String excep = "G4IonisationParameters - data file: " 
                     + name + " not found";
      G4Exception(excep);
    }

    // The file is organized into two columns:
    // 1st column is the energy
    // 2nd column is the corresponding value
    // The file terminates with the pattern: -1   -1
    //                                       -2   

    e   = new G4DataVector();
    for(size_t j=0; j<length; j++) {
      a[j] = new G4DataVector();
    } 
    for(size_t k=0; k<length; k++) {
      G4int id = Z*20 + k;
      p[k] = new G4ShellEMDataSet(id, interpolation, 1., 1.);
    }

    do {
      file >> energy;
      if (energy == -2) break;
      if (energy >  -1) e->push_back(energy);

      G4double x;
      for (size_t j=0; j<length; j++) {
          file >> x;    
          if(energy > -1) a[j]->push_back(x);
      }

      // fill map
      if(energy < 0) {

        G4int id;
        for(G4int k=0; k<length; k++) {
            set = new G4EMDataSet(k, e, a[k], interpolation, 1., 1.);
            p[k]->AddComponent(set);
            G4int id = Z*20 + k;
        } 

	// clear vectors
        e   = new G4DataVector();
        for(size_t j2=0; j2<length; j2++) {
          a[j2] = new G4DataVector();
        } 
      }
    } while (energy > -2);

    file.close();

    for(size_t kk=0; kk<length; kk++) {
      G4int id = Z*20 + kk;
      param[id] = p[kk];
    }

    // remove last vectors because they are not used
    delete e;
    for(size_t j1=0; j1<length; j1++) {
      delete  a[j1];
    } 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eIonisationParameters::PrintData() const
{
  G4cout << G4endl;
  G4cout << "===== G4eIonisationParameters =====" << G4endl;
  G4cout << G4endl;

  size_t nZ = activeZ.size();
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;

  for (G4int i=0; i<nZ; i++) {
    G4int Z = (G4int)activeZ[i];      

    for (G4int j=0; j<length; j++) {
    
      G4int index = Z*20 + j;

      pos = param.find(index);
      if (pos!= param.end()) {
        G4VEMDataSet* dataSet = pos->second;
        size_t nShells = dataSet->NumberOfComponents();

        for (G4int k=0; k<nShells; k++) {

          G4cout << "===== Z= " << Z << " shell= " << k 
                 << " parameter[" << j << "]  =====" 
                 << G4endl;
          const G4VEMDataSet* comp = dataSet->GetComponent(k);
          comp->PrintData();
	}
      }
    }
  }
  G4cout << "====================================" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

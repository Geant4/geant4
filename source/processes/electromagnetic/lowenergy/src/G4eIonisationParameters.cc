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
// $Id: G4eIonisationParameters.cc,v 1.9 2001-10-24 22:02:18 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
// 12.09.01 V.Ivanchenko    Add param and interpolation of parameters  
// 04.10.01 V.Ivanchenko    Add BindingEnergy method  
//
// -------------------------------------------------------------------

#include "G4eIonisationParameters.hh"
#include "G4VEMDataSet.hh"
#include "G4ShellEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4Material.hh"
#include "G4DataVector.hh"
#include "g4std/fstream"
#include "g4std/strstream"


G4eIonisationParameters:: G4eIonisationParameters(G4int minZ, G4int maxZ)
  : zMin(minZ), zMax(maxZ),
  length(16)
{
  LoadData();
}


G4eIonisationParameters::~G4eIonisationParameters()
{ 
  // Reset the map of data sets: remove the data sets from the map 
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::iterator pos;

  for (pos = param.begin(); pos != param.end(); ++pos)
    {
      G4VEMDataSet* dataSet = (*pos).second;
      //      param.erase(pos);
      delete dataSet;
      dataSet = 0;
    }

  activeZ.clear();
}


G4double G4eIonisationParameters::Parameter(G4int Z, G4int shellIndex, 
                                                     G4int parameterIndex, 
                                                     G4double e) const
{
  G4double value = 0.;
  G4int id = Z*20 + parameterIndex;
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;

  pos = param.find(id);
  if (pos!= param.end()) {
    G4VEMDataSet* dataSet = (*pos).second;
    G4int nShells = dataSet->NumberOfComponents();

    if(shellIndex < nShells) { 
      const G4VEMDataSet* comp = dataSet->GetComponent(shellIndex);
      const G4DataVector ener = comp->GetEnergies(0);
      G4double ee = G4std::max(ener.front(),G4std::min(ener.back(),e));
      value = comp->FindValue(ee);

    } else {
      G4cout << "WARNING: G4IonisationParameters::FindParameters "
             << "has no parameters for shell= " << shellIndex 
             << "; Z= " << Z
             << G4endl;
    }
  } else {
    G4cout << "WARNING: G4IonisationParameters::FindValue "
           << "did not find ID = "
           << shellIndex << G4endl;
  }

  return value;
}


void G4eIonisationParameters::LoadData()
{
  // ---------------------------------------
  // Please document what are the parameters 
  // ---------------------------------------

  // define active elements

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
     G4Exception("G4CrossSectionHandler: no MaterialTable found)");

  G4int nMaterials = G4Material::GetNumberOfMaterials();
  
  for (G4int m=0; m<nMaterials; m++) {

    const G4Material* material= (*materialTable)[m];        
    const G4ElementVector* elementVector = material->GetElementVector();
    const size_t nElements = material->GetNumberOfElements();
      
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
  G4DataVector  q;
  G4VEMDataSet* set;
  G4std::vector<G4DataVector*> a;
  G4std::vector<G4ShellEMDataSet*> p;
  a.resize(length);
  p.resize(length);
  q.resize(length);
  
  size_t nZ = activeZ.size();
  
  for (size_t i=0; i<nZ; i++) {
    
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
      G4VDataSetAlgorithm* interpolation  = new G4SemiLogInterpolation();
      p[k] = new G4CompositeEMDataSet(id, interpolation, 1., 1.);
    }

    G4int shell = 0;

    do {
      file >> energy;
      if (energy == -2) break;
      if (energy >  -1) e->push_back(energy);

      size_t j;
      for (j=0; j<length; j++) {
        file >> q[j];
      }    
      if(q[14] == 0.)   q[14] = DBL_MAX;
      if(q[15] < q[14]) q[15] = q[14];
      if(q[15] == 0.)   q[15] = DBL_MAX;

      for (j=0; j<length; j++) {
        a[j]->push_back(q[j]);
      }

      // fill map
      if(energy < 0) {

        for(size_t k=0; k<length; k++) {
            G4int id = Z*20 + k;
	    G4VDataSetAlgorithm* interpolation  = new G4SemiLogInterpolation();
            set = new G4EMDataSet(id, e, a[k], interpolation, 1., 1.);
            p[k]->AddComponent(set);
        } 

	// clear vectors
        e   = new G4DataVector();
        for(size_t j2=0; j2<length; j2++) {
          a[j2] = new G4DataVector();
        } 
        shell++;
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
  a.clear();
  p.clear();
}


void G4eIonisationParameters::PrintData() const
{
  G4cout << G4endl;
  G4cout << "===== G4eIonisationParameters =====" << G4endl;
  G4cout << G4endl;

  size_t nZ = activeZ.size();
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;

  for (size_t i=0; i<nZ; i++) {
    G4int Z = (G4int)activeZ[i];      

    for (size_t j=0; j<length; j++) {
    
      G4int index = Z*20 + j;

      pos = param.find(index);
      if (pos!= param.end()) {
        G4VEMDataSet* dataSet = (*pos).second;
        size_t nShells = dataSet->NumberOfComponents();

        for (size_t k=0; k<nShells; k++) {

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



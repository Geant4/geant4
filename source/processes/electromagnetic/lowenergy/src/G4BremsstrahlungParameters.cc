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
// $Id: G4BremsstrahlungParameters.cc,v 1.10 2001-11-29 19:01:36 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//         V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
// 12.09.01 V.Ivanchenko    Add activeZ and paramA
// 25.09.01 V.Ivanchenko    Add parameter C and change interface to B
// 29.11.01 V.Ivanchenko    Update parametrisation
//
// -------------------------------------------------------------------

#include "G4BremsstrahlungParameters.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4LogLogInterpolation.hh"
#include "G4Material.hh"
#include "g4std/fstream"
#include "g4std/strstream"


G4BremsstrahlungParameters:: G4BremsstrahlungParameters(G4int minZ, G4int maxZ)
  : zMin(minZ), 
    zMax(maxZ),
    length(16)
{
  LoadData();
}


G4BremsstrahlungParameters::~G4BremsstrahlungParameters()
{ 
  // Reset the map of data sets: remove the data sets from the map 
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::iterator pos;

  for (pos = param.begin(); pos != param.end(); pos++)
    {
      G4VEMDataSet* dataSet = (*pos).second;
      delete dataSet;
    }

  activeZ.clear();
  paramC.clear();
}


G4double G4BremsstrahlungParameters::Parameter(G4int parameterIndex, 
                                                 G4int Z, 
                                                 G4double energy) const
{
  G4double value = 0.;
  G4int id = Z*20 + parameterIndex;
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;

  pos = param.find(id);
  if (pos!= param.end()) {

    G4VEMDataSet* dataSet = (*pos).second;
    const G4DataVector ener = dataSet->GetEnergies(0);
    G4double ee = G4std::max(ener.front(),G4std::min(ener.back(),energy));
    value = dataSet->FindValue(ee);
    
  } else {
    G4cout << "WARNING: G4BremsstrahlungParameters::FindValue "
           << "did not find ID = "
           << id << G4endl;
  }

  return value;
}

void G4BremsstrahlungParameters::LoadData()
{
  // Build the complete string identifying the file with the data set

  // define active elements

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
     G4Exception("G4CrossSectionHandler: no MaterialTable found)");

  G4int nMaterials = G4Material::GetNumberOfMaterials();

  G4double x = 1.e-9;
  for (G4int mm=0; mm<100; mm++) {
    paramC.push_back(x);
  }
  
  for (G4int m=0; m<nMaterials; m++) {

    const G4Material* material= (*materialTable)[m];        
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4int nElements = material->GetNumberOfElements();
      
    for (G4int iEl=0; iEl<nElements; iEl++) {
      G4Element* element = (*elementVector)[iEl];
      G4double Z = element->GetZ();
      G4int iz = (G4int)Z;
      if(iz < 100)
            paramC[iz] = 0.217635e-33*(material->GetTotNbOfElectPerVolume());
      if (!(activeZ.contains(Z))) {
	 activeZ.push_back(Z);
      }   
    }
  }

  // Read parameters
   
  char* path = getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep = G4String("G4BremsstrahlungParameters - G4LEDATA") 
                     + G4String("environment variable not set");
      G4Exception(excep);
    }

  G4String pathString_a(path);
  G4String name_a = pathString_a + "/brem/br-sp.dat";  
  G4std::ifstream file_a(name_a);
  G4std::filebuf* lsdp_a = file_a.rdbuf();
  
  if (! (lsdp_a->is_open()) ) {
     G4String excep = G4String("G4BremsstrahlungParameters: cannot open file ")
                    + name_a;
     G4Exception(excep);
  }  

  // The file is organized into two columns:
  // 1st column is the energy
  // 2nd column is the corresponding value
  // The file terminates with the pattern: -1   -1
  //                                       -2   -2

  G4DataVector* energies;
  G4DataVector* data;
  G4double ener = 0.0;
  G4double sum = 0.0;
  size_t nZ  = activeZ.size();
  energies   = new G4DataVector(); 
  data       = new G4DataVector(); 
  G4int z    = 0;

  G4bool used = false;
  G4std::vector<G4DataVector*> a;
  for (size_t j=0; j<length; j++) {
    G4DataVector* aa = new G4DataVector();
    a.push_back(aa);
  } 
  G4DataVector e;
  e.clear();

  do {
    file_a >> ener >> sum;

    // End of file
    if (ener == -2) {
      break;

      // End of next element
    } else if (ener == -1) {

      z++;
      G4double Z = (G4double)z;
    
	// fill map if Z is used
      if (activeZ.contains(Z)) {

 	for (size_t k=0; k<length; k++) {

	    G4int id = z*20 + k;
	    G4VDataSetAlgorithm* inter  = new G4LogLogInterpolation();
	    G4DataVector* eVector = new G4DataVector;
	    size_t eSize = e.size();
	    for (size_t s=0; s<eSize; s++) {
	       eVector->push_back(e[s]);
	    }
	    G4VEMDataSet* set = new G4EMDataSet(id,eVector,a[k],inter,1.,1.);
    	    param[id] = set;
	}
        used = true;
        a.clear();
        for (size_t j=0; j<length; j++) {
            G4DataVector* aa = new G4DataVector();
            a.push_back(aa);
        } 
      }
      if(!used) {
        for (size_t j=0; j<length; j++) {
          a[j]->clear();
          used = false; 
        }
      }
      e.clear();
  
    } else {

      if(ener > 1000.) ener = 1000.;
      e.push_back(ener);
      a[length-1]->push_back(sum);

      for (size_t j=0; j<length-1; j++) {
	G4double qRead;
	file_a >> qRead;
	/*
        if(ener == 1000.) {
          G4double x = 0.1*((G4double)j);
          if(j == 0) x = 0.01;
          if(j == 10) x = 0.95;
          if(j == 11) x = 0.97;
          if(j == 12) x = 0.99;
          if(j == 13) x = 0.995;
          if(j == 14) x = 1.0;
          qRead = 1. - x + 0.75*x*x;
	}
	*/
	a[j]->push_back(qRead);
      }    

    }
  } while (ener != -2);
  
  file_a.close();

}


G4double G4BremsstrahlungParameters::ParameterC(G4int id) const
{
  G4int n = paramC.size();
  if (id < 0 || id >= n) {
    G4String ex = "G4BremsstrahlungParameters::ParameterC - wrong id=" + id; 
    G4Exception(ex);
  }

  return paramC[id];
}


void G4BremsstrahlungParameters::PrintData() const
{
  
  G4cout << G4endl;
  G4cout << "===== G4BremsstrahlungParameters =====" << G4endl;
  G4cout << G4endl;
  G4cout << "===== Parameters =====" << G4endl;
  G4cout << G4endl;

  size_t nZ = activeZ.size();
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;

  for (size_t j=0; j<nZ; j++) {
    G4int Z = (G4int)activeZ[j];   

    for (size_t i=0; i<length; i++) {

      pos = param.find(Z*20 + i);
      if (pos!= param.end()) {

        G4cout << "===== Z= " << Z  
                 << " parameter[" << i << "]  =====" 
                 << G4endl;
        G4VEMDataSet* dataSet = (*pos).second;
        dataSet->PrintData();
      }
    }
  }

  G4cout << "==========================================" << G4endl;
}


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
// $Id: G4BremsstrahlungParameters.cc,v 1.2 2001-09-14 08:17:13 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
// 13.09.01 V.Ivanchenko  Add paramA map
//
// -------------------------------------------------------------------

#include "G4BremsstrahlungParameters.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4LogLogInterpolation.hh"
#include "G4MaterialTable.hh"
#include "g4std/fstream"
#include "g4std/strstream"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BremsstrahlungParameters:: G4BremsstrahlungParameters(G4int minZ, G4int maxZ)
  : zMin(minZ), zMax(maxZ)
{
  LoadData();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BremsstrahlungParameters::~G4BremsstrahlungParameters()
{ 
  // Reset the map of data sets: remove the data sets from the map 
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::iterator pos;

  for (pos = paramA.begin(); pos != paramA.end(); pos++)
    {
      G4VEMDataSet* dataSet = pos->second;
      paramA.erase(pos);
      delete dataSet;
    }

  activeZ.clear();
  delete interpolation;
  paramB0.clear();
  paramB1.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BremsstrahlungParameters::ParameterA(G4int Z, G4double energy) const
{
  G4double value = 0.;
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;
  pos = paramA.find(Z);
  if (pos!= paramA.end()) {
    G4VEMDataSet* dataSet = pos->second;
    value = dataSet->FindValue(energy);
  } else {
    G4String ex = "G4BremsstrahlungParameters::ParameterA - wrong Z=" + Z; 
    G4Exception(ex);
  }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BremsstrahlungParameters::ParameterB(G4int Z, G4int id) const
{
  if (Z < (zMin-1) || Z > (zMax-1)) {
    G4String ex = "G4BremsstrahlungParameters::ParameterB - wrong Z=" + Z; 
    G4Exception(ex);
  }

  G4double value;

  if (id == 1) 
    {
      value = paramB0[Z-1];
    }
  else
    {
      value = paramB1[Z-1];
    }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BremsstrahlungParameters::LoadData()
{
  // Build the complete string identifying the file with the data set

  if(!interpolation) interpolation = new G4LogLogInterpolation();

  // define active elements

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
     G4Exception("G4CrossSectionHandler: no MaterialTable found)");

  G4int nMaterials = materialTable->length();
  
  for (G4int m=0; m<nMaterials; m++) {

    const G4Material* material= (*materialTable)[m];        
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4int nElements = material->GetNumberOfElements();
      
    for (G4int iEl=0; iEl<nElements; iEl++) {
      G4Element* element = (*elementVector)[iEl];
      G4double Z = element->GetZ();
      if (!(activeZ.contains(Z)) && Z > 0 && Z < zMax) {
	 activeZ.push_back(Z);
      }
    }
  }

  // Read parameters A
   
  char* path = getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep = G4String("G4BremsstrahlungParameters - G4LEDATA") +
                     + G4String("environment variable not set");
      G4Exception(excep);
    }

  G4String pathString_a(path);
  G4String name_a = pathString_a + "/brem/br-co-a.dat";  
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
  G4double e = 0;
  G4double a = 0;
  size_t nZ  = activeZ.size();
  energies   = new G4DataVector(); 
  data       = new G4DataVector(); 
  G4int z    = 0;

  do {
    file_a >> e >> a;

    // End of the next element
    if (e == -1) {

      z++;
      G4bool used = false;
      for (size_t i=0; i<nZ; i++) {
    
	// fill map
        if(z == (G4int)activeZ[i]) {
          G4VEMDataSet* dataSet = 
                        new G4EMDataSet(z, energies, data, interpolation);
          paramA[z] = dataSet;
          used = true;
	}
      }
      if(!used) {
        energies->clear();
        data->clear();
      }
  
      if (z > 98) break;

      energies   = new G4DataVector(); 
      data       = new G4DataVector(); 
      
    } else {

      energies->push_back(e);
      data->push_back(a);
    }
  } while (e != -2);
  
  file_a.close();

  // Read parameters B
  
  char nameChar[100] = {""};
  G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
  
  G4String fileName = "brem/br-co-b";
  ost << fileName  << ".dat";
  
  G4String name(nameChar);
    
  G4String pathString(path);
  G4String dirFile = pathString + "/" + name;
  G4std::ifstream file(dirFile);
  G4std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
	{
	  G4String excep = "G4BremsstrahlungParameters - data file: " + dirFile + " not found";
	  G4Exception(excep);
	}

  G4int nParam = 2;
  a = 0;
  G4int k = 1;

  do
    {
      file >> a;
      // The file is organized into two columns:
      // 1st column is the energy
      // 2nd column is the corresponding value
      // The file terminates with the pattern: -1   -1
      //                                       -2   -2

      if (k%nParam != 0)
	{	
	  paramB0.push_back(a);
	  k++;
	}
      else if (k%nParam == 0)
	{
	  paramB1.push_back(a);
	  k = 1;
	}
    } while (a != -2); // end of file
  
  zMax = paramB0.size();
  file.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BremsstrahlungParameters::PrintData() const
{
  
  G4cout << G4endl;
  G4cout << "===== G4BremsstrahlungParameters =====" << G4endl;
  G4cout << G4endl;
  G4cout << "===== Parameter A =====" << G4endl;
  G4cout << G4endl;

  size_t nZ = activeZ.size();
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;

  for (G4int j=0; j<nZ; j++) {
    G4int Z = (G4int)activeZ[j];   
    pos = paramA.find(Z);
    if (pos!= paramA.end()) {
      G4cout << "===== Z= " << Z << " =====" << G4endl;
      G4VEMDataSet* dataSet = pos->second;
      dataSet->PrintData();
    }
  }

  G4cout << G4endl;
  G4cout << "===== Parameter B =====" << G4endl;
  G4cout << G4endl;

  for (G4int i=zMin-1; i<zMax-1; i++)
    {
      G4double b0 = paramB0[i];
      G4double b1 = paramB1[i];
      G4cout << "Z = " << i << " - Parameter B: "
	     << b0
	     << ",  "
	     << b1
	     << G4endl; 
    }
  G4cout << "==========================================" << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

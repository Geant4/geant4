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
// $Id: G4BremsstrahlungParameters.cc,v 1.6 2001-10-24 22:02:18 pia Exp $
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
    zMax(maxZ)
{
  LoadData();
}


G4BremsstrahlungParameters::~G4BremsstrahlungParameters()
{ 
  // Reset the map of data sets: remove the data sets from the map 
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::iterator pos;

  for (pos = paramA.begin(); pos != paramA.end(); pos++)
    {
      G4VEMDataSet* dataSet = (*pos).second;
      paramA.erase(pos);
      delete dataSet;
    }

  activeZ.clear();
  paramB0.clear();
  paramB1.clear();
  paramC.clear();
}


G4double G4BremsstrahlungParameters::ParameterA(G4int Z, G4double energy) const
{
  G4double value = 0.;
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;
  pos = paramA.find(Z);
  if (pos!= paramA.end()) {
    G4VEMDataSet* dataSet = (*pos).second;
    value = dataSet->FindValue(energy);
  } else {
    G4String ex = "G4BremsstrahlungParameters::ParameterA - wrong Z=" + Z; 
    G4Exception(ex);
  }
  return value;
}


G4double G4BremsstrahlungParameters::ParameterB(G4int Z, G4double energy) const
{
  if (Z < zMin || Z > zMax) {
    G4String ex = "G4BremsstrahlungParameters::ParameterB - wrong Z=" + Z; 
    G4Exception(ex);
  }

  G4double b0  = paramB0[Z-1];
  G4double b1  = paramB1[Z-1];
  G4double b   = b0 + b1*log10(energy);
  if(b < 0. || energy < 10.*eV) b = 0.;

  return b;
}


void G4BremsstrahlungParameters::LoadData()
{
  // Build the complete string identifying the file with the data set

  // define active elements

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
     G4Exception("G4CrossSectionHandler: no MaterialTable found)");

  G4int nMaterials = G4Material::GetNumberOfMaterials();
  
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
    G4double x = 1.e-9;
    paramC.push_back(x);
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
        if (z == (G4int)activeZ[i]) {
	  G4VDataSetAlgorithm* interpolation  = new G4LogLogInterpolation();
          G4VEMDataSet* dataSet = 
                    new G4EMDataSet(z, energies, data, interpolation, 1., 1.);
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

  //  G4int nParam = 2;
  a = 0;
  e = 0;

  do {
      file >> e >> a;
      // The file is organized into two columns:
      // 1st column is the energy
      // 2nd column is the corresponding value
      // The file terminates with the pattern: -1   -1
      //                                       -2   -2

      if (e == -1 || e == -2) break;
      else {
	  paramB0.push_back(e);
	  paramB1.push_back(a);
	}
    } while (e != -2); // end of file
  
  zMax = paramB0.size();
  file.close();
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
  G4cout << "===== Parameter A =====" << G4endl;
  G4cout << G4endl;

  size_t nZ = activeZ.size();
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;

  for (size_t j=0; j<nZ; j++) {
    G4int Z = (G4int)activeZ[j];   
    pos = paramA.find(Z);
    if (pos!= paramA.end()) {
      G4cout << "===== Z= " << Z << " =====" << G4endl;
      G4VEMDataSet* dataSet = (*pos).second;
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


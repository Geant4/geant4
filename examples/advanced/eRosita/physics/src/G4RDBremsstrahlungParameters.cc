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
//
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
// 18.11.02 V.Ivanchenko    Fix problem of load
// 21.02.03 V.Ivanchenko    Number of parameters is defined in the constructor
// 28.02.03 V.Ivanchenko    Filename is defined in the constructor
//
// -------------------------------------------------------------------

#include "G4RDBremsstrahlungParameters.hh"
#include "G4RDVEMDataSet.hh"
#include "G4RDEMDataSet.hh"
#include "G4RDLogLogInterpolation.hh"
#include "G4Material.hh"
#include <fstream>


G4RDBremsstrahlungParameters:: G4RDBremsstrahlungParameters(const G4String& name,
    size_t num, G4int minZ, G4int maxZ)
  : zMin(minZ),
    zMax(maxZ),
    length(num)
{
  LoadData(name);
}


G4RDBremsstrahlungParameters::~G4RDBremsstrahlungParameters()
{
  // Reset the map of data sets: remove the data sets from the map
  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> >::iterator pos;

  for (pos = param.begin(); pos != param.end(); ++pos)
    {
      G4RDVEMDataSet* dataSet = (*pos).second;
      delete dataSet;
    }

  activeZ.clear();
  paramC.clear();
}


G4double G4RDBremsstrahlungParameters::Parameter(G4int parameterIndex,
                                                 G4int Z,
                                                 G4double energy) const
{
  G4double value = 0.;
  G4int id = Z*length + parameterIndex;
  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> >::const_iterator pos;

  pos = param.find(id);
  if (pos!= param.end()) {

    G4RDVEMDataSet* dataSet = (*pos).second;
    const G4DataVector ener = dataSet->GetEnergies(0);
    G4double ee = std::max(ener.front(),std::min(ener.back(),energy));
    value = dataSet->FindValue(ee);

  } else {
    G4cout << "WARNING: G4RDBremsstrahlungParameters::FindValue "
           << "did not find ID = "
           << id << G4endl;
  }

  return value;
}

void G4RDBremsstrahlungParameters::LoadData(const G4String& name)
{
  // Build the complete string identifying the file with the data set

  // define active elements

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
     G4Exception("G4RDBremsstrahlungParameters::LoadData()",
                 "DataNotFound", FatalException, "No MaterialTable found!");

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

  char* path = std::getenv("G4LEDATA");
  if (path == 0)
    {
      G4String excep("G4LEDATA environment variable not set!");
      G4Exception("G4RDBremsstrahlungParameters::LoadData()",
                  "InvalidSetup", FatalException, excep);
    }

  G4String pathString_a(path);
  G4String name_a = pathString_a + name;
  std::ifstream file_a(name_a);
  std::filebuf* lsdp_a = file_a.rdbuf();

  if (! (lsdp_a->is_open()) )
    {
      G4String stringConversion2("Cannot open file ");
      G4String excep = stringConversion2 + name_a;
      G4Exception("G4RDBremsstrahlungParameters::LoadData()",
                  "FileNotFound", FatalException, excep);
  }

  // The file is organized into two columns:
  // 1st column is the energy
  // 2nd column is the corresponding value
  // The file terminates with the pattern: -1   -1
  //                                       -2   -2

  G4double ener = 0.0;
  G4double sum = 0.0;
  G4int z    = 0;

  std::vector<G4DataVector*> a;
  for (size_t j=0; j<length; j++) {
    G4DataVector* aa = new G4DataVector();
    a.push_back(aa);
  }
  G4DataVector e;
  e.clear();

  do {
    file_a >> ener >> sum;

    // End of file
    if (ener == (G4double)(-2)) {
      break;

      // End of next element
    } else if (ener == (G4double)(-1)) {

      z++;
      G4double Z = (G4double)z;

	// fill map if Z is used
      if (activeZ.contains(Z)) {

 	for (size_t k=0; k<length; k++) {

	    G4int id = z*length + k;
	    G4RDVDataSetAlgorithm* inter  = new G4RDLogLogInterpolation();
	    G4DataVector* eVector = new G4DataVector;
	    size_t eSize = e.size();
	    for (size_t s=0; s<eSize; s++) {
	       eVector->push_back(e[s]);
	    }
	    G4RDVEMDataSet* set = new G4RDEMDataSet(id,eVector,a[k],inter,1.,1.);
    	    param[id] = set;
	}
        a.clear();
        for (size_t j=0; j<length; j++) {
            G4DataVector* aa = new G4DataVector();
            a.push_back(aa);
        }
      } else {
        for (size_t j=0; j<length; j++) {
          a[j]->clear();
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
	a[j]->push_back(qRead);
      }

    }
  } while (ener != (G4double)(-2));

  file_a.close();

}


G4double G4RDBremsstrahlungParameters::ParameterC(G4int id) const
{
  G4int n = paramC.size();
  if (id < 0 || id >= n)
    {
      G4String stringConversion1("Wrong id = ");
      G4String stringConversion2(id);
      G4String ex = stringConversion1 + stringConversion2;
      G4Exception("G4RDBremsstrahlungParameters::ParameterC()",
                  "InvalidSetup", FatalException, ex);
    }

  return paramC[id];
}


void G4RDBremsstrahlungParameters::PrintData() const
{

  G4cout << G4endl;
  G4cout << "===== G4RDBremsstrahlungParameters =====" << G4endl;
  G4cout << G4endl;
  G4cout << "===== Parameters =====" << G4endl;
  G4cout << G4endl;

  size_t nZ = activeZ.size();
  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> >::const_iterator pos;

  for (size_t j=0; j<nZ; j++) {
    G4int Z = (G4int)activeZ[j];

    for (size_t i=0; i<length; i++) {

      pos = param.find(Z*length + i);
      if (pos!= param.end()) {

        G4cout << "===== Z= " << Z
                 << " parameter[" << i << "]  ====="
                 << G4endl;
        G4RDVEMDataSet* dataSet = (*pos).second;
        dataSet->PrintData();
      }
    }
  }

  G4cout << "==========================================" << G4endl;
}


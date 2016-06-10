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
// $Id: G4BremsstrahlungParameters.cc 66241 2012-12-13 18:34:42Z gunter $
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
// 03.12.10 V.Ivanchenko    Fixed memory leak in LoadData
//
// -------------------------------------------------------------------

#include <fstream>

#include "G4BremsstrahlungParameters.hh"
#include "G4PhysicalConstants.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4LogLogInterpolation.hh"
#include "G4Material.hh"

G4BremsstrahlungParameters:: G4BremsstrahlungParameters(const G4String& name,
    size_t num, G4int minZ, G4int maxZ)
  : zMin(minZ),
    zMax(maxZ),
    length(num)
{
  LoadData(name);
}


G4BremsstrahlungParameters::~G4BremsstrahlungParameters()
{
  // Reset the map of data sets: remove the data sets from the map
  std::map<G4int,G4VEMDataSet*,std::less<G4int> >::iterator pos;

  for (pos = param.begin(); pos != param.end(); ++pos)
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
  G4int id = Z*length + parameterIndex;
  std::map<G4int,G4VEMDataSet*,std::less<G4int> >::const_iterator pos;

  pos = param.find(id);
  if (pos!= param.end()) {

    G4VEMDataSet* dataSet = (*pos).second;
    const G4DataVector ener = dataSet->GetEnergies(0);
    G4double ee = std::max(ener.front(),std::min(ener.back(),energy));
    value = dataSet->FindValue(ee);

  } else {
    G4cout << "WARNING: G4BremsstrahlungParameters::FindValue "
           << "did not find ID = "
           << id << G4endl;
  }

  return value;
}

void G4BremsstrahlungParameters::LoadData(const G4String& name)
{
  const G4double mConst = 
    classic_electr_radius*electron_Compton_length*electron_Compton_length*4.0*pi;

  // Build the complete string identifying the file with the data set
  // define active elements

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
     G4Exception("G4BremsstrahlungParameters::LoadData",
		    "em1001",FatalException,"Unable to find MaterialTable");

  G4int nMaterials = G4Material::GetNumberOfMaterials();

  G4double x = 1.e-9;
  for (G4int mmLocal=0; mmLocal<100; mmLocal++) {
    paramC.push_back(x);
  }

  for (G4int mLocal=0; mLocal<nMaterials; mLocal++) {

    const G4Material* material= (*materialTable)[mLocal];
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4int nElements = material->GetNumberOfElements();

    for (G4int iEl=0; iEl<nElements; iEl++) {
      G4Element* element = (*elementVector)[iEl];
      G4double Z = element->GetZ();
      G4int iz = (G4int)Z;
      if(iz < 100) {
	paramC[iz] = mConst*material->GetTotNbOfElectPerVolume();
	//paramC[iz] = 0.217635e-33*(material->GetTotNbOfElectPerVolume());
      }
      if (!(activeZ.contains(Z))) {
	 activeZ.push_back(Z);
      }
    }
  }

  // Read parameters

  char* path = getenv("G4LEDATA");
  if (path == 0)
    {
      G4Exception("G4BremsstrahlungParameters::LoadData",
		    "em0006",FatalException,"G4LEDATA environment variable not set");
      return;
    }

  G4String pathString_a(path);
  G4String name_a = pathString_a + name;
  std::ifstream file_a(name_a);
  std::filebuf* lsdp_a = file_a.rdbuf();

  if (! (lsdp_a->is_open()) )
    {
      G4String stringConversion2("G4BremsstrahlungParameters::LoadData");
      G4String excep = stringConversion2 + name_a;
      G4Exception("G4BremsstrahlungParameters::LoadData",
		    "em0003",FatalException,excep);
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
  a.resize(length);
 
  G4DataVector e;
  e.clear();

  G4bool isReady = false;

  do {
    file_a >> ener >> sum;

    // End of file
    if (ener == (G4double)(-2)) {
      break;

      // End of next element
    } else if (ener == (G4double)(-1)) {

      ++z;
      G4double Z = (G4double)z;

      // fill map if Z is used
      if (activeZ.contains(Z)) {

 	for (size_t k=0; k<length; ++k) {

	  G4int id = z*length + k;
	  G4VDataSetAlgorithm* inter  = new G4LogLogInterpolation();
	  G4DataVector* eVector = new G4DataVector;
	  size_t eSize = e.size();
	  for (size_t sLocal=0; sLocal<eSize; sLocal++) {
	    eVector->push_back(e[sLocal]);
	  }
	  G4VEMDataSet* set = new G4EMDataSet(id,eVector,a[k],inter,1.,1.);
	  param[id] = set;
	}
      } else {
        for (size_t j=0; j<length; j++) {
	  delete a[j];
	}
      }
      isReady = false;

    } else {

      if(!isReady) {
        isReady = true;
	e.clear();
        for (size_t j=0; j<length; ++j) {
          a[j] = new G4DataVector();
        }
      }

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


G4double G4BremsstrahlungParameters::ParameterC(G4int id) const
{
  G4int n = paramC.size();
  if (id < 0 || id >= n)
    {
      G4String stringConversion2(id);
      G4String ex = "Wrong id " + stringConversion2;
      G4Exception("G4BremsstrahlungParameters::ParameterC",
		    "em1002",FatalException,ex);

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
  std::map<G4int,G4VEMDataSet*,std::less<G4int> >::const_iterator pos;

  for (size_t j=0; j<nZ; j++) {
    G4int Z = (G4int)activeZ[j];

    for (size_t i=0; i<length; i++) {

      pos = param.find(Z*length + i);
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


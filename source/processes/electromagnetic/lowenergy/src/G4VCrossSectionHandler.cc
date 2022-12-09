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
//
// History:
// -----------
// 1  Aug 2001   MGP        Created
// 09 Oct 2001   VI         Add FindValue with 3 parameters 
//                          + NumberOfComponents
// 19 Jul 2002   VI         Create composite data set for material
// 21 Jan 2003   VI         Cut per region
//
// 15 Jul 2009   Nicolas A. Karakatsanis
//
//                           - LoadNonLogData method was created to load only the non-logarithmic data from G4EMLOW
//                             dataset. It is essentially performing the data loading operations as in the past.
//
//                           - LoadData method was revised in order to calculate the logarithmic values of the data
//                             It retrieves the data values from the G4EMLOW data files but, then, calculates the
//                             respective log values and loads them to seperate data structures.
//                             The EM data sets, initialized this way, contain both non-log and log values.
//                             These initialized data sets can enhance the computing performance of data interpolation
//                             operations
//
//                           - BuildMeanFreePathForMaterials method was also revised in order to calculate the 
//                             logarithmic values of the loaded data. 
//                             It generates the data values and, then, calculates the respective log values which 
//                             later load to seperate data structures.
//                             The EM data sets, initialized this way, contain both non-log and log values.
//                             These initialized data sets can enhance the computing performance of data interpolation
//                             operations
//                             
//                           - LoadShellData method was revised in order to eliminate the presence of a potential
//                             memory leak originally identified by Riccardo Capra.
//                             Riccardo Capra Original Comment
//                             Riccardo Capra <capra@ge.infn.it>: PLEASE CHECK THE FOLLOWING PIECE OF CODE
//                             "energies" AND "data" G4DataVector ARE ALLOCATED, FILLED IN AND NEVER USED OR
//                             DELETED. WHATSMORE LOADING FILE OPERATIONS WERE DONE BY G4ShellEMDataSet
//                             EVEN BEFORE THE CHANGES I DID ON THIS FILE. SO THE FOLLOWING CODE IN MY
//                             OPINION SHOULD BE USELESS AND SHOULD PRODUCE A MEMORY LEAK.
//
//
// -------------------------------------------------------------------

#include "G4VCrossSectionHandler.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4ShellEMDataSet.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "Randomize.hh"
#include <map>
#include <vector>
#include <fstream>
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VCrossSectionHandler::G4VCrossSectionHandler()
{
  crossSections = 0;
  interpolation = 0;
  Initialise();
  ActiveElements();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VCrossSectionHandler::G4VCrossSectionHandler(G4VDataSetAlgorithm* algorithm,
					       G4double minE,
					       G4double maxE,
					       G4int bins,
					       G4double unitE,
					       G4double unitData,
					       G4int minZ, 
					       G4int maxZ)
  : interpolation(algorithm), eMin(minE), eMax(maxE), 
    unit1(unitE), unit2(unitData), zMin(minZ), zMax(maxZ), 
    nBins(bins)
{
  crossSections = nullptr;
  ActiveElements();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VCrossSectionHandler::~G4VCrossSectionHandler()
{
  delete interpolation;
  interpolation = nullptr;

  for (auto & pos : dataMap)
    {
      G4VEMDataSet* dataSet = pos.second;
      delete dataSet;
    }

  if (crossSections != nullptr)
    {
      std::size_t n = crossSections->size();
      for (std::size_t i=0; i<n; i++)
	{
	  delete (*crossSections)[i];
	}
      delete crossSections;
      crossSections = nullptr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VCrossSectionHandler::Initialise(G4VDataSetAlgorithm* algorithm,
					G4double minE, G4double maxE, 
					G4int numberOfBins,
					G4double unitE, G4double unitData,
					G4int minZ, G4int maxZ)
{
  if (algorithm != nullptr) 
    {
      delete interpolation;
      interpolation = algorithm;
    }
  else
    {
      delete interpolation;
      interpolation = CreateInterpolation();
    }

  eMin = minE;
  eMax = maxE;
  nBins = numberOfBins;
  unit1 = unitE;
  unit2 = unitData;
  zMin = minZ;
  zMax = maxZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VCrossSectionHandler::PrintData() const
{
  for (auto & pos : dataMap) 
    {
      G4int z = pos.first;
      G4VEMDataSet* dataSet = pos.second;     
      G4cout << "---- Data set for Z = "
	     << z
	     << G4endl;
      dataSet->PrintData();
      G4cout << "--------------------------------------------------" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VCrossSectionHandler::LoadData(const G4String& fileName)
{
  std::size_t nZ = activeZ.size();
  for (std::size_t i=0; i<nZ; ++i)
    {
      G4int Z = G4int(activeZ[i]);

      // Build the complete string identifying the file with the data set      
      const char* path = G4FindDataDir("G4LEDATA");
      if (!path)
	{ 
          G4Exception("G4VCrossSectionHandler::LoadData",
		    "em0006",FatalException,"G4LEDATA environment variable not set");
	  return;
	}
      
      std::ostringstream ost;
      ost << path << '/' << fileName << Z << ".dat";
      std::ifstream file(ost.str().c_str());
      std::filebuf* lsdp = file.rdbuf();
       
      if (! (lsdp->is_open()) )
	{
	  G4String excep = "data file: " + ost.str() + " not found";
          G4Exception("G4VCrossSectionHandler::LoadData",
		    "em0003",FatalException,excep);
	}
      G4double a = 0;
      G4int k = 0;
      G4int nColumns = 2;

      G4DataVector* orig_reg_energies = new G4DataVector;
      G4DataVector* orig_reg_data = new G4DataVector;
      G4DataVector* log_reg_energies = new G4DataVector;
      G4DataVector* log_reg_data = new G4DataVector;

      do
	{
	  file >> a;

          if (a==0.) a=1e-300;

	  // The file is organized into four columns:
	  // 1st column contains the values of energy
	  // 2nd column contains the corresponding data value
	  // The file terminates with the pattern: -1   -1
	  //                                       -2   -2
          //
	  if (a != -1 && a != -2)
	    {
	      if (k%nColumns == 0)
                {
		 orig_reg_energies->push_back(a*unit1);
                 log_reg_energies->push_back(std::log10(a)+std::log10(unit1));
                }
	      else if (k%nColumns == 1)
                {
		 orig_reg_data->push_back(a*unit2);
                 log_reg_data->push_back(std::log10(a)+std::log10(unit2));
                }
              k++;
	    }
	} 
      while (a != -2); // End of File
      
      file.close();
      G4VDataSetAlgorithm* algo = interpolation->Clone();
      G4VEMDataSet* dataSet = new G4EMDataSet(Z,orig_reg_energies,orig_reg_data,
					      log_reg_energies,log_reg_data,algo);
      dataMap[Z] = dataSet;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VCrossSectionHandler::LoadNonLogData(const G4String& fileName)
{
  std::size_t nZ = activeZ.size();
  for (std::size_t i=0; i<nZ; ++i)
    {
      G4int Z = G4int(activeZ[i]);

      // Build the complete string identifying the file with the data set      
      const char* path = G4FindDataDir("G4LEDATA");
      if (!path)
	{ 
          G4Exception("G4VCrossSectionHandler::LoadNonLogData",
		    "em0006",FatalException,"G4LEDATA environment variable not set");
	  return;
	}
      
      std::ostringstream ost;
      ost << path << '/' << fileName << Z << ".dat";
      std::ifstream file(ost.str().c_str());
      std::filebuf* lsdp = file.rdbuf();
       
      if (! (lsdp->is_open()) )
	{
	  G4String excep = "data file: " + ost.str() + " not found";
          G4Exception("G4VCrossSectionHandler::LoadNonLogData",
		    "em0003",FatalException,excep);
	}
      G4double a = 0;
      G4int k = 0;
      G4int nColumns = 2;

      G4DataVector* orig_reg_energies = new G4DataVector;
      G4DataVector* orig_reg_data = new G4DataVector;

      do
	{
	  file >> a;

	  // The file is organized into four columns:
	  // 1st column contains the values of energy
	  // 2nd column contains the corresponding data value
	  // The file terminates with the pattern: -1   -1
	  //                                       -2   -2
          //
	  if (a != -1 && a != -2)
	    {
	      if (k%nColumns == 0)
                {
		 orig_reg_energies->push_back(a*unit1);
                }
	      else if (k%nColumns == 1)
                {
		 orig_reg_data->push_back(a*unit2);
                }
              k++;
	    }
	} 
      while (a != -2); // End of File
      
      file.close();
      G4VDataSetAlgorithm* algo = interpolation->Clone();

      G4VEMDataSet* dataSet = new G4EMDataSet(Z,orig_reg_energies,orig_reg_data,algo);
      dataMap[Z] = dataSet;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VCrossSectionHandler::LoadShellData(const G4String& fileName)
{
  std::size_t nZ = activeZ.size();
  for (std::size_t i=0; i<nZ; ++i)
    {
      G4int Z = G4int(activeZ[i]);
      
      G4VDataSetAlgorithm* algo = interpolation->Clone();
      G4VEMDataSet* dataSet = new G4ShellEMDataSet(Z, algo);
      dataSet->LoadData(fileName);      
      dataMap[Z] = dataSet;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VCrossSectionHandler::Clear()
{
  // Reset the map of data sets: remove the data sets from the map 
  if(! dataMap.empty())
    {
      for (auto & pos : dataMap)
	{
	  G4VEMDataSet* dataSet = pos.second;
	  delete dataSet;
	  dataSet = nullptr;
	  G4int i = pos.first;
	  dataMap[i] = nullptr;
	}
	dataMap.clear();
    }
  activeZ.clear();
  ActiveElements();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VCrossSectionHandler::FindValue(G4int Z, G4double energy) const
{
  G4double value = 0.;
 
  auto pos = dataMap.find(Z);
  if (pos!= dataMap.end())
    {
      G4VEMDataSet* dataSet = (*pos).second;
      value = dataSet->FindValue(energy);
    }
  else
    {
      G4cout << "WARNING: G4VCrossSectionHandler::FindValue did not find Z = "
	     << Z << G4endl;
    }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VCrossSectionHandler::FindValue(G4int Z, G4double energy, 
                                           G4int shellIndex) const
{
  G4double value = 0.;
  auto pos = dataMap.find(Z);
  if (pos!= dataMap.cend())
    {
      G4VEMDataSet* dataSet = (*pos).second;
      if (shellIndex >= 0) 
	{
	  G4int nComponents = (G4int)dataSet->NumberOfComponents();
	  if(shellIndex < nComponents)    
	    value = dataSet->GetComponent(shellIndex)->FindValue(energy);
	  else 
	    {
	      G4cout << "WARNING: G4VCrossSectionHandler::FindValue did not find"
		     << " shellIndex= " << shellIndex
		     << " for  Z= "
		     << Z << G4endl;
	    }
	} else {
	  value = dataSet->FindValue(energy);
	}
    }
  else
    {
      G4cout << "WARNING: G4VCrossSectionHandler::FindValue did not find Z = "
	     << Z << G4endl;
    }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VCrossSectionHandler::ValueForMaterial(const G4Material* material,
						  G4double energy) const
{
  G4double value = 0.;
  const G4ElementVector* elementVector = material->GetElementVector();
  const G4double* nAtomsPerVolume = material->GetVecNbOfAtomsPerVolume();
  std::size_t nElements = material->GetNumberOfElements();

  for (std::size_t i=0 ; i<nElements ; ++i)
    {
      G4int Z = (G4int) (*elementVector)[i]->GetZ();
      G4double elementValue = FindValue(Z,energy);
      G4double nAtomsVol = nAtomsPerVolume[i];
      value += nAtomsVol * elementValue;
    }

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEMDataSet* G4VCrossSectionHandler::BuildMeanFreePathForMaterials(const G4DataVector* energyCuts)
{
  // Builds a CompositeDataSet containing the mean free path for each material
  // in the material table
  G4DataVector energyVector;
  G4double dBin = std::log10(eMax/eMin) / nBins;

  for (G4int i=0; i<nBins+1; i++)
    {
      energyVector.push_back(std::pow(10., std::log10(eMin)+i*dBin));
    }

  // Factory method to build cross sections in derived classes,
  // related to the type of physics process

  if (crossSections != nullptr)
    {  // Reset the list of cross sections
      if (! crossSections->empty())
	{
	  for (auto mat=crossSections->begin(); mat != crossSections->end(); ++mat)
	    {
	      G4VEMDataSet* set = *mat;
	      delete set;
	      set = nullptr;
	    }
	  crossSections->clear();
	  delete crossSections;
	  crossSections = nullptr;
	}
    }

  crossSections = BuildCrossSectionsForMaterials(energyVector,energyCuts);

  if (crossSections == nullptr)
    {
      G4Exception("G4VCrossSectionHandler::BuildMeanFreePathForMaterials",
		    "em1010",FatalException,"crossSections = 0");
      return 0;
    }

  G4VDataSetAlgorithm* algo = CreateInterpolation();
  G4VEMDataSet* materialSet = new G4CompositeEMDataSet(algo);

  G4DataVector* energies;
  G4DataVector* data;
  G4DataVector* log_energies;
  G4DataVector* log_data;
  
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

  for (G4int mLocal=0; mLocal<numOfCouples; ++mLocal)
    {
      energies = new G4DataVector;
      data = new G4DataVector;
      log_energies = new G4DataVector;
      log_data = new G4DataVector;
      for (G4int bin=0; bin<nBins; bin++)
	{
	  G4double energy = energyVector[bin];
	  energies->push_back(energy);
          log_energies->push_back(std::log10(energy));
	  G4VEMDataSet* matCrossSet = (*crossSections)[mLocal];
	  G4double materialCrossSection = 0.0;
          G4int nElm = (G4int)matCrossSet->NumberOfComponents();
          for(G4int j=0; j<nElm; ++j) {
            materialCrossSection += matCrossSet->GetComponent(j)->FindValue(energy);
	  }

	  if (materialCrossSection > 0.)
	    {
	      data->push_back(1./materialCrossSection);
              log_data->push_back(std::log10(1./materialCrossSection));
	    }
	  else
	    {
	      data->push_back(DBL_MAX);
              log_data->push_back(std::log10(DBL_MAX));
	    }
	}
      G4VDataSetAlgorithm* algoLocal = CreateInterpolation();
      G4VEMDataSet* dataSet = new G4EMDataSet(mLocal,energies,data,log_energies,log_data,algoLocal,1.,1.);
      materialSet->AddComponent(dataSet);
    }
  return materialSet;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4VCrossSectionHandler::SelectRandomAtom(const G4MaterialCutsCouple* couple,
                                                     G4double e) const
{
  // Select randomly an element within the material, according to the weight
  // determined by the cross sections in the data set
  const G4Material* material = couple->GetMaterial();
  G4int nElements = (G4int)material->GetNumberOfElements();

  // Special case: the material consists of one element
  if (nElements == 1)
    {
      G4int Z = (G4int) material->GetZ();
      return Z;
    }

  // Composite material
  const G4ElementVector* elementVector = material->GetElementVector();
  std::size_t materialIndex = couple->GetIndex();

  G4VEMDataSet* materialSet = (*crossSections)[materialIndex];
  G4double materialCrossSection0 = 0.0;
  G4DataVector cross;
  cross.clear();
  for ( G4int i=0; i < nElements; i++ )
    {
      G4double cr = materialSet->GetComponent(i)->FindValue(e);
      materialCrossSection0 += cr;
      cross.push_back(materialCrossSection0);
    }

  G4double random = G4UniformRand() * materialCrossSection0;
  for (G4int k=0 ; k < nElements ; k++ )
    {
      if (random <= cross[k]) return (G4int) (*elementVector)[k]->GetZ();
    }
  // It should never get here
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Element* G4VCrossSectionHandler::SelectRandomElement(const G4MaterialCutsCouple* couple,
							     G4double e) const
{
  // Select randomly an element within the material, according to the weight determined
  // by the cross sections in the data set
  const G4Material* material = couple->GetMaterial();
  G4Element* nullElement = 0;
  G4int nElements = (G4int)material->GetNumberOfElements();
  const G4ElementVector* elementVector = material->GetElementVector();

  // Special case: the material consists of one element
  if (nElements == 1)
    {
      return (*elementVector)[0];
    }
  else
    {
      // Composite material

      std::size_t materialIndex = couple->GetIndex();

      G4VEMDataSet* materialSet = (*crossSections)[materialIndex];
      G4double materialCrossSection0 = 0.0;
      G4DataVector cross;
      cross.clear();
      for (G4int i=0; i<nElements; ++i)
        {
          G4double cr = materialSet->GetComponent(i)->FindValue(e);
          materialCrossSection0 += cr;
          cross.push_back(materialCrossSection0);
        }

      G4double random = G4UniformRand() * materialCrossSection0;

      for (G4int k=0 ; k < nElements ; ++k )
        {
          if (random <= cross[k]) return (*elementVector)[k];
        }
      // It should never end up here
      G4cout << "G4VCrossSectionHandler::SelectRandomElement - no element found" << G4endl;
      return nullElement;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4VCrossSectionHandler::SelectRandomShell(G4int Z, G4double e) const
{
  // Select randomly a shell, according to the weight determined by the cross sections
  // in the data set
  // Note for later improvement: it would be useful to add a cache mechanism for already
  // used shells to improve performance
  G4int shell = 0;

  G4double totCrossSection = FindValue(Z,e);
  G4double random = G4UniformRand() * totCrossSection;
  G4double partialSum = 0.;

  G4VEMDataSet* dataSet = nullptr;
  auto pos = dataMap.find(Z);
  if (pos != dataMap.end()) 
    dataSet = (*pos).second;
  else
    {
      G4Exception("G4VCrossSectionHandler::SelectRandomShell",
		    "em1011",FatalException,"unable to load the dataSet");
      return 0;
    }

  G4int nShells = (G4int)dataSet->NumberOfComponents();
  for (G4int i=0; i<nShells; ++i)
    {
      const G4VEMDataSet* shellDataSet = dataSet->GetComponent(i);
      if (shellDataSet != nullptr)
	{
	  G4double value = shellDataSet->FindValue(e);
	  partialSum += value;
	  if (random <= partialSum) return i;
	}
    }
  // It should never get here
  return shell;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VCrossSectionHandler::ActiveElements()
{
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == nullptr)
      G4Exception("G4VCrossSectionHandler::ActiveElements",
		    "em1001",FatalException,"no MaterialTable found");

  std::size_t nMaterials = G4Material::GetNumberOfMaterials();

  for (std::size_t mLocal2=0; mLocal2<nMaterials; ++mLocal2)
    {
      const G4Material* material= (*materialTable)[mLocal2];
      const G4ElementVector* elementVector = material->GetElementVector();
      const std::size_t nElements = material->GetNumberOfElements();

      for (std::size_t iEl=0; iEl<nElements; ++iEl)
	{
	  G4double Z = (*elementVector)[iEl]->GetZ();
	  if (!(activeZ.contains(Z)) && Z >= zMin && Z <= zMax)
	    {
	      activeZ.push_back(Z);
	    }
	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VDataSetAlgorithm* G4VCrossSectionHandler::CreateInterpolation()
{
  G4VDataSetAlgorithm* algorithm = new G4LogLogInterpolation;
  return algorithm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4VCrossSectionHandler::NumberOfComponents(G4int Z) const
{
  G4int n = 0;

  auto pos = dataMap.find(Z);
  if (pos!= dataMap.end())
    {
      G4VEMDataSet* dataSet = (*pos).second;
      n = (G4int)dataSet->NumberOfComponents();
    }
  else
    {
      G4cout << "WARNING: G4VCrossSectionHandler::NumberOfComponents did not "
             << "find Z = "
             << Z << G4endl;
    }
  return n;
}



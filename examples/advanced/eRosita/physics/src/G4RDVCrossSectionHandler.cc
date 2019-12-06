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
// -------------------------------------------------------------------

#include "G4RDVCrossSectionHandler.hh"
#include "G4RDVDataSetAlgorithm.hh"
#include "G4RDLogLogInterpolation.hh"
#include "G4RDVEMDataSet.hh"
#include "G4RDEMDataSet.hh"
#include "G4RDCompositeEMDataSet.hh"
#include "G4RDShellEMDataSet.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "Randomize.hh"
#include <map>
#include <vector>
#include <fstream>
#include <sstream>


G4RDVCrossSectionHandler::G4RDVCrossSectionHandler()
{
  crossSections = 0;
  interpolation = 0;
  Initialise();
  ActiveElements();
}


G4RDVCrossSectionHandler::G4RDVCrossSectionHandler(G4RDVDataSetAlgorithm* algorithm,
					       G4double minE,
					       G4double maxE,
					       G4int bins,
					       G4double unitE,
					       G4double unitData,
					       G4int minZ, 
					       G4int maxZ)
  : interpolation(algorithm), eMin(minE), eMax(maxE), nBins(bins),
    unit1(unitE), unit2(unitData), zMin(minZ), zMax(maxZ)
{
  crossSections = 0;
  ActiveElements();
}

G4RDVCrossSectionHandler::~G4RDVCrossSectionHandler()
{
  delete interpolation;
  interpolation = 0;
  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> >::iterator pos;

  for (pos = dataMap.begin(); pos != dataMap.end(); ++pos)
    {
      // The following is a workaround for STL ObjectSpace implementation, 
      // which does not support the standard and does not accept 
      // the syntax pos->second
      // G4RDVEMDataSet* dataSet = pos->second;
      G4RDVEMDataSet* dataSet = (*pos).second;
      delete dataSet;
    }

  if (crossSections != 0)
    {
      size_t n = crossSections->size();
      for (size_t i=0; i<n; i++)
	{
	  delete (*crossSections)[i];
	}
      delete crossSections;
      crossSections = 0;
    }
}

void G4RDVCrossSectionHandler::Initialise(G4RDVDataSetAlgorithm* algorithm,
					G4double minE, G4double maxE, 
					G4int numberOfBins,
					G4double unitE, G4double unitData,
					G4int minZ, G4int maxZ)
{
  if (algorithm != 0) 
    {
      delete interpolation;
      interpolation = algorithm;
    }
  else
    {
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

void G4RDVCrossSectionHandler::PrintData() const
{
  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> >::const_iterator pos;

  for (pos = dataMap.begin(); pos != dataMap.end(); pos++)
    {
      // The following is a workaround for STL ObjectSpace implementation, 
      // which does not support the standard and does not accept 
      // the syntax pos->first or pos->second
      // G4int z = pos->first;
      // G4RDVEMDataSet* dataSet = pos->second;
      G4int z = (*pos).first;
      G4RDVEMDataSet* dataSet = (*pos).second;     
      G4cout << "---- Data set for Z = "
	     << z
	     << G4endl;
      dataSet->PrintData();
      G4cout << "--------------------------------------------------" << G4endl;
    }
}

void G4RDVCrossSectionHandler::LoadData(const G4String& fileName)
{
  size_t nZ = activeZ.size();
  for (size_t i=0; i<nZ; i++)
    {
      G4int Z = (G4int) activeZ[i];

      // Build the complete string identifying the file with the data set
      
      char* path = std::getenv("G4LEDATA");
      if (!path)
	{ 
	  G4String excep = "G4LEDATA environment variable not set!";
	  G4Exception("G4RDVCrossSectionHandler::LoadData()",
                      "InvalidSetup", FatalException, excep);
	}
      
      std::ostringstream ost;
      ost << path << '/' << fileName << Z << ".dat";
      std::ifstream file(ost.str().c_str());
      std::filebuf* lsdp = file.rdbuf();
      
      if (! (lsdp->is_open()) )
	{
	  G4String excep = "Data file: " + ost.str() + " not found!";
	  G4Exception("G4RDVCrossSectionHandler::LoadData()",
                      "DataNotFound", FatalException, excep);
	}
      G4double a = 0;
      G4int k = 1;
      G4DataVector* energies = new G4DataVector;
      G4DataVector* data = new G4DataVector;
      do
	{
	  file >> a;
	  G4int nColumns = 2;
	  // The file is organized into two columns:
	  // 1st column is the energy
	  // 2nd column is the corresponding value
	  // The file terminates with the pattern: -1   -1
	  //                                       -2   -2
	  if (a == -1 || a == -2)
	    {
	    }
	  else
	    {
	      if (k%nColumns != 0)
		{	
		  G4double e = a * unit1;
		  energies->push_back(e);
		  k++;
		}
	      else if (k%nColumns == 0)
		{
		  G4double value = a * unit2;
		  data->push_back(value);
		  k = 1;
		}
	    }
	} while (a != -2); // end of file
      
      file.close();
      G4RDVDataSetAlgorithm* algo = interpolation->Clone();
      G4RDVEMDataSet* dataSet = new G4RDEMDataSet(Z,energies,data,algo);
      dataMap[Z] = dataSet;
    }
}

void G4RDVCrossSectionHandler::LoadShellData(const G4String& fileName)
{
  size_t nZ = activeZ.size();
  for (size_t i=0; i<nZ; i++)
    {
      G4int Z = (G4int) activeZ[i];
      
      // Riccardo Capra <capra@ge.infn.it>: PLEASE CHECK THE FOLLOWING PIECE OF CODE
      // "energies" AND "data" G4DataVector ARE ALLOCATED, FILLED IN AND NEVER USED OR
      // DELETED. WHATSMORE LOADING FILE OPERATIONS WERE DONE BY G4RDShellEMDataSet
      // EVEN BEFORE THE CHANGES I DID ON THIS FILE. SO THE FOLLOWING CODE IN MY
      // OPINION SHOULD BE USELESS AND SHOULD PRODUCE A MEMORY LEAK. 

      // Build the complete string identifying the file with the data set
      
      char* path = std::getenv("G4LEDATA");
      if (!path)
	{ 
	  G4String excep = "G4LEDATA environment variable not set!";
	  G4Exception("G4RDVCrossSectionHandler::LoadShellData()",
                      "InvalidSetup", FatalException, excep);
	}
      
      std::ostringstream ost;

      ost << path << '/' << fileName << Z << ".dat";
      
      std::ifstream file(ost.str().c_str());
      std::filebuf* lsdp = file.rdbuf();
      
      if (! (lsdp->is_open()) )
	{
	  G4String excep = "Data file: " + ost.str() + " not found!";
	  G4Exception("G4RDVCrossSectionHandler::LoadShellData()",
                      "DataNotFound", FatalException, excep);
	}
      G4double a = 0;
      G4int k = 1;
      G4DataVector* energies = new G4DataVector;
      G4DataVector* data = new G4DataVector;
      do
	{
	  file >> a;
	  G4int nColumns = 2;
	  // The file is organized into two columns:
	  // 1st column is the energy
	  // 2nd column is the corresponding value
	  // The file terminates with the pattern: -1   -1
	  //                                       -2   -2
	  if (a == -1 || a == -2)
	    {
	    }
	  else
	    {
	      if (k%nColumns != 0)
		{	
		  G4double e = a * unit1;
		  energies->push_back(e);
		  k++;
		}
	      else if (k%nColumns == 0)
		{
		  G4double value = a * unit2;
		  data->push_back(value);
		  k = 1;
		}
	    }
	} while (a != -2); // end of file
      
      file.close();
      
      // Riccardo Capra <capra@ge.infn.it>: END OF CODE THAT IN MY OPINION SHOULD BE
      // REMOVED.
      
      G4RDVDataSetAlgorithm* algo = interpolation->Clone();
      G4RDVEMDataSet* dataSet = new G4RDShellEMDataSet(Z, algo);
      dataSet->LoadData(fileName);
      dataMap[Z] = dataSet;
    }
}

void G4RDVCrossSectionHandler::Clear()
{
  // Reset the map of data sets: remove the data sets from the map 
  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> >::iterator pos;

  if(! dataMap.empty())
    {
        for (pos = dataMap.begin(); pos != dataMap.end(); ++pos)
	{
	  // The following is a workaround for STL ObjectSpace implementation, 
	  // which does not support the standard and does not accept
	  // the syntax pos->first or pos->second
	  // G4RDVEMDataSet* dataSet = pos->second;
	  G4RDVEMDataSet* dataSet = (*pos).second;
	  delete dataSet;
	  dataSet = 0;
	  G4int i = (*pos).first;
	  dataMap[i] = 0;
	}
	dataMap.clear();
    }

  activeZ.clear();
  ActiveElements();
}

G4double G4RDVCrossSectionHandler::FindValue(G4int Z, G4double energy) const
{
  G4double value = 0.;
  
  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> >::const_iterator pos;
  pos = dataMap.find(Z);
  if (pos!= dataMap.end())
    {
      // The following is a workaround for STL ObjectSpace implementation, 
      // which does not support the standard and does not accept 
      // the syntax pos->first or pos->second
      // G4RDVEMDataSet* dataSet = pos->second;
      G4RDVEMDataSet* dataSet = (*pos).second;
      value = dataSet->FindValue(energy);
    }
  else
    {
      G4cout << "WARNING: G4RDVCrossSectionHandler::FindValue did not find Z = "
	     << Z << G4endl;
    }
  return value;
}

G4double G4RDVCrossSectionHandler::FindValue(G4int Z, G4double energy, 
                                           G4int shellIndex) const
{
  G4double value = 0.;

  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> >::const_iterator pos;
  pos = dataMap.find(Z);
  if (pos!= dataMap.end())
    {
      // The following is a workaround for STL ObjectSpace implementation, 
      // which does not support the standard and does not accept 
      // the syntax pos->first or pos->second
      // G4RDVEMDataSet* dataSet = pos->second;
      G4RDVEMDataSet* dataSet = (*pos).second;
      if (shellIndex >= 0) 
	{
	  G4int nComponents = dataSet->NumberOfComponents();
	  if(shellIndex < nComponents)    
	    // - MGP - Why doesn't it use G4RDVEMDataSet::FindValue directly?
	    value = dataSet->GetComponent(shellIndex)->FindValue(energy);
	  else 
	    {
	      G4cout << "WARNING: G4RDVCrossSectionHandler::FindValue did not find"
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
      G4cout << "WARNING: G4RDVCrossSectionHandler::FindValue did not find Z = "
	     << Z << G4endl;
    }
  return value;
}


G4double G4RDVCrossSectionHandler::ValueForMaterial(const G4Material* material,
						  G4double energy) const
{
  G4double value = 0.;

  const G4ElementVector* elementVector = material->GetElementVector();
  const G4double* nAtomsPerVolume = material->GetVecNbOfAtomsPerVolume();
  G4int nElements = material->GetNumberOfElements();

  for (G4int i=0 ; i<nElements ; i++)
    {
      G4int Z = (G4int) (*elementVector)[i]->GetZ();
      G4double elementValue = FindValue(Z,energy);
      G4double nAtomsVol = nAtomsPerVolume[i];
      value += nAtomsVol * elementValue;
    }

  return value;
}


G4RDVEMDataSet* G4RDVCrossSectionHandler::BuildMeanFreePathForMaterials(const G4DataVector* energyCuts)
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

  if (crossSections != 0)
    {  // Reset the list of cross sections
      std::vector<G4RDVEMDataSet*>::iterator mat;
      if (! crossSections->empty())
	{
	  for (mat = crossSections->begin(); mat!= crossSections->end(); ++mat)
	    {
	      G4RDVEMDataSet* set = *mat;
	      delete set;
	      set = 0;
	    }
	  crossSections->clear();
	  delete crossSections;
	  crossSections = 0;
	}
    }

  crossSections = BuildCrossSectionsForMaterials(energyVector,energyCuts);

  if (crossSections == 0)
    G4Exception("G4RDVCrossSectionHandler::BuildMeanFreePathForMaterials()",
                "InvalidCondition", FatalException, "CrossSections = 0!");

  G4RDVDataSetAlgorithm* algo = CreateInterpolation();
  G4RDVEMDataSet* materialSet = new G4RDCompositeEMDataSet(algo);

  G4DataVector* energies;
  G4DataVector* data;
  
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();


  for (size_t m=0; m<numOfCouples; m++)
    {
      energies = new G4DataVector;
      data = new G4DataVector;
      for (G4int bin=0; bin<nBins; bin++)
	{
	  G4double energy = energyVector[bin];
	  energies->push_back(energy);
	  G4RDVEMDataSet* matCrossSet = (*crossSections)[m];
	  G4double materialCrossSection = 0.0;
          G4int nElm = matCrossSet->NumberOfComponents();
          for(G4int j=0; j<nElm; j++) {
            materialCrossSection += matCrossSet->GetComponent(j)->FindValue(energy);
	  }

	  if (materialCrossSection > 0.)
	    {
	      data->push_back(1./materialCrossSection);
	    }
	  else
	    {
	      data->push_back(DBL_MAX);
	    }
	}
      G4RDVDataSetAlgorithm* algo = CreateInterpolation();
      G4RDVEMDataSet* dataSet = new G4RDEMDataSet(m,energies,data,algo,1.,1.);
      materialSet->AddComponent(dataSet);
    }

  return materialSet;
}

G4int G4RDVCrossSectionHandler::SelectRandomAtom(const G4MaterialCutsCouple* couple,
                                                     G4double e) const
{
  // Select randomly an element within the material, according to the weight
  // determined by the cross sections in the data set

  const G4Material* material = couple->GetMaterial();
  G4int nElements = material->GetNumberOfElements();

  // Special case: the material consists of one element
  if (nElements == 1)
    {
      G4int Z = (G4int) material->GetZ();
      return Z;
    }

  // Composite material

  const G4ElementVector* elementVector = material->GetElementVector();
  size_t materialIndex = couple->GetIndex();

  G4RDVEMDataSet* materialSet = (*crossSections)[materialIndex];
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

const G4Element* G4RDVCrossSectionHandler::SelectRandomElement(const G4MaterialCutsCouple* couple,
							     G4double e) const
{
  // Select randomly an element within the material, according to the weight determined
  // by the cross sections in the data set

  const G4Material* material = couple->GetMaterial();
  G4Element* nullElement = 0;
  G4int nElements = material->GetNumberOfElements();
  const G4ElementVector* elementVector = material->GetElementVector();

  // Special case: the material consists of one element
  if (nElements == 1)
    {
      G4Element* element = (*elementVector)[0];
      return element;
    }
  else
    {
      // Composite material

      size_t materialIndex = couple->GetIndex();

      G4RDVEMDataSet* materialSet = (*crossSections)[materialIndex];
      G4double materialCrossSection0 = 0.0;
      G4DataVector cross;
      cross.clear();
      for (G4int i=0; i<nElements; i++)
        {
          G4double cr = materialSet->GetComponent(i)->FindValue(e);
          materialCrossSection0 += cr;
          cross.push_back(materialCrossSection0);
        }

      G4double random = G4UniformRand() * materialCrossSection0;

      for (G4int k=0 ; k < nElements ; k++ )
        {
          if (random <= cross[k]) return (*elementVector)[k];
        }
      // It should never end up here
      G4cout << "G4RDVCrossSectionHandler::SelectRandomElement - no element found" << G4endl;
      return nullElement;
    }
}

G4int G4RDVCrossSectionHandler::SelectRandomShell(G4int Z, G4double e) const
{
  // Select randomly a shell, according to the weight determined by the cross sections
  // in the data set

  // Note for later improvement: it would be useful to add a cache mechanism for already
  // used shells to improve performance

  G4int shell = 0;

  G4double totCrossSection = FindValue(Z,e);
  G4double random = G4UniformRand() * totCrossSection;
  G4double partialSum = 0.;

  G4RDVEMDataSet* dataSet = 0;
  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> >::const_iterator pos;
  pos = dataMap.find(Z);
  // The following is a workaround for STL ObjectSpace implementation,
  // which does not support the standard and does not accept
  // the syntax pos->first or pos->second
  // if (pos != dataMap.end()) dataSet = pos->second;
  if (pos != dataMap.end()) dataSet = (*pos).second;

  size_t nShells = dataSet->NumberOfComponents();
  for (size_t i=0; i<nShells; i++)
    {
      const G4RDVEMDataSet* shellDataSet = dataSet->GetComponent(i);
      if (shellDataSet != 0)
	{
	  G4double value = shellDataSet->FindValue(e);
	  partialSum += value;
	  if (random <= partialSum) return i;
	}
    }
  // It should never get here
  return shell;
}

void G4RDVCrossSectionHandler::ActiveElements()
{
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
    G4Exception("G4RDVCrossSectionHandler::ActiveElements",
                "InvalidSetup", FatalException, "No MaterialTable found!");

  G4int nMaterials = G4Material::GetNumberOfMaterials();

  for (G4int m=0; m<nMaterials; m++)
    {
      const G4Material* material= (*materialTable)[m];
      const G4ElementVector* elementVector = material->GetElementVector();
      const G4int nElements = material->GetNumberOfElements();

      for (G4int iEl=0; iEl<nElements; iEl++)
	{
	  G4Element* element = (*elementVector)[iEl];
	  G4double Z = element->GetZ();
	  if (!(activeZ.contains(Z)) && Z >= zMin && Z <= zMax)
	    {
	      activeZ.push_back(Z);
	    }
	}
    }
}

G4RDVDataSetAlgorithm* G4RDVCrossSectionHandler::CreateInterpolation()
{
  G4RDVDataSetAlgorithm* algorithm = new G4RDLogLogInterpolation;
  return algorithm;
}

G4int G4RDVCrossSectionHandler::NumberOfComponents(G4int Z) const
{
  G4int n = 0;

  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> >::const_iterator pos;
  pos = dataMap.find(Z);
  if (pos!= dataMap.end())
    {
      G4RDVEMDataSet* dataSet = (*pos).second;
      n = dataSet->NumberOfComponents();
    }
  else
    {
      G4cout << "WARNING: G4RDVCrossSectionHandler::NumberOfComponents did not "
             << "find Z = "
             << Z << G4endl;
    }
  return n;
}



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
// $Id: G4VCrossSectionHandler.cc,v 1.6 2001-10-09 11:23:28 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 1 Aug 2001   MGP        Created
// 09.10.01  V.Ivanchenko  Add FindValue with 3 parameters 
//                          + NumberOfComponents
//
// -------------------------------------------------------------------

#include "G4VCrossSectionHandler.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4ShellEMDataSet.hh"
#include "G4MaterialTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "Randomize.hh" 
#include "g4std/map"
#include "g4std/vector"
#include "g4std/fstream"
#include "g4std/strstream"


G4VCrossSectionHandler::G4VCrossSectionHandler()
{
  crossSections = 0;
  interpolation = 0;
  Initialise();
  ActiveElements();
}


G4VCrossSectionHandler::G4VCrossSectionHandler(G4VDataSetAlgorithm* algorithm,
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

G4VCrossSectionHandler::~G4VCrossSectionHandler()
{
  delete interpolation;
  interpolation = 0;
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::iterator pos;

  for (pos = dataMap.begin(); pos != dataMap.end(); ++pos)
    {
      // The following is a workaround for STL ObjectSpace implementation, 
      // which does not support the standard and does not accept 
      // the syntax pos->second
      // G4VEMDataSet* dataSet = pos->second;
      G4VEMDataSet* dataSet = (*pos).second;
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

void G4VCrossSectionHandler::Initialise(G4VDataSetAlgorithm* algorithm,
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

void G4VCrossSectionHandler::PrintData() const
{
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;

  for (pos = dataMap.begin(); pos != dataMap.end(); pos++)
    {
      // The following is a workaround for STL ObjectSpace implementation, 
      // which does not support the standard and does not accept 
      // the syntax pos->first or pos->second
      // G4int z = pos->first;
      // G4VEMDataSet* dataSet = pos->second;
      G4int z = (*pos).first;
      G4VEMDataSet* dataSet = (*pos).second;     
      G4cout << "---- Data set for Z = "
	     << z
	     << G4endl;
      dataSet->PrintData();
      G4cout << "--------------------------------------------------" << G4endl;
    }
}

void G4VCrossSectionHandler::LoadData(const G4String& fileName)
{
  size_t nZ = activeZ.size();
  for (size_t i=0; i<nZ; i++)
    {
      G4int Z = (G4int) activeZ[i];

      // Build the complete string identifying the file with the data set
      
      char nameChar[100] = {""};
      G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
      
      ost << fileName << Z << ".dat";
      
      G4String name(nameChar);
      
      char* path = getenv("G4LEDATA");
      if (!path)
	{ 
	  G4String excep = "G4VCrossSectionHandler - G4LEDATA environment variable not set";
	  G4Exception(excep);
	}
      
      G4String pathString(path);
      G4String dirFile = pathString + "/" + name;
      G4std::ifstream file(dirFile);
      G4std::filebuf* lsdp = file.rdbuf();
      
      if (! (lsdp->is_open()) )
	{
	  G4String excep = "G4VCrossSectionHandler - data file: " + dirFile + " not found";
	  G4Exception(excep);
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
      G4VDataSetAlgorithm* algo = interpolation->Clone();
      G4VEMDataSet* dataSet = new G4EMDataSet(Z,energies,data,algo);
      dataMap[Z] = dataSet;
    }
}

void G4VCrossSectionHandler::LoadShellData(const G4String& fileName)
{
  size_t nZ = activeZ.size();
  for (size_t i=0; i<nZ; i++)
    {
      G4int Z = (G4int) activeZ[i];

      // Build the complete string identifying the file with the data set
      
      char nameChar[100] = {""};
      G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
      
      ost << fileName << Z << ".dat";
      
      G4String name(nameChar);
      
      char* path = getenv("G4LEDATA");
      if (!path)
	{ 
	  G4String excep = "G4VCrossSectionHandler - G4LEDATA environment variable not set";
	  G4Exception(excep);
	}
      
      G4String pathString(path);
      G4String dirFile = pathString + "/" + name;
      G4std::ifstream file(dirFile);
      G4std::filebuf* lsdp = file.rdbuf();
      
      if (! (lsdp->is_open()) )
	{
	  G4String excep = "G4VCrossSectionHandler - data file: " + dirFile + " not found";
	  G4Exception(excep);
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
      G4VDataSetAlgorithm* algo = interpolation->Clone();
      G4VEMDataSet* dataSet = new G4ShellEMDataSet(Z,fileName,algo);
      dataMap[Z] = dataSet;
    }
}

void G4VCrossSectionHandler::Clear()
{
  // Reset the map of data sets: remove the data sets from the map 
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::iterator pos;

  if(! dataMap.empty())
    {
        for (pos = dataMap.begin(); pos != dataMap.end(); ++pos)
	{
	  // The following is a workaround for STL ObjectSpace implementation, 
	  // which does not support the standard and does not accept 
	  // the syntax pos->first or pos->second
	  // G4VEMDataSet* dataSet = pos->second;
	  G4VEMDataSet* dataSet = (*pos).second;
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

G4double G4VCrossSectionHandler::FindValue(G4int Z, G4double energy) const
{
  G4double value = 0.;
  
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;
  pos = dataMap.find(Z);
  if (pos!= dataMap.end())
    {
      // The following is a workaround for STL ObjectSpace implementation, 
      // which does not support the standard and does not accept 
      // the syntax pos->first or pos->second
      // G4VEMDataSet* dataSet = pos->second;
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

G4double G4VCrossSectionHandler::FindValue(G4int Z, G4double energy, 
                                           G4int shellIndex) const
{
  G4double value = 0.;
  
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;
  pos = dataMap.find(Z);
  if (pos!= dataMap.end())
    {
      // The following is a workaround for STL ObjectSpace implementation, 
      // which does not support the standard and does not accept 
      // the syntax pos->first or pos->second
      // G4VEMDataSet* dataSet = pos->second;
      G4VEMDataSet* dataSet = (*pos).second;
      if(shellIndex >= 0) {
        if(shellIndex < dataSet->NumberOfComponents())
          value = dataSet->GetComponent(shellIndex)->FindValue(energy);
        else {
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


G4double G4VCrossSectionHandler::ValueForMaterial(const G4Material* material, 
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


G4VEMDataSet* G4VCrossSectionHandler::BuildMeanFreePathForMaterials(
                                      const G4DataVector* energyCuts) 
{
  // Builds a CompositeDataSet containing the mean free path for each material
  // in the material table

  G4DataVector energyVector;
  G4double dBin = log10(eMax/eMin) / nBins;

  for (G4int i=0; i<nBins+1; i++) 
    {
      energyVector.push_back(pow(10., log10(eMin)+i*dBin));
    }

  // Factory method to build cross sections in derived classes, 
  // related to the type of physics process 

  if (crossSections != 0)
    {  // Reset the list of cross sections
      G4std::vector<G4VEMDataSet*>::iterator mat;
      if (! crossSections->empty())
	{
	  for (mat = crossSections->begin(); mat!= crossSections->end(); ++mat)
	    {
	      G4VEMDataSet* set = *mat;
	      delete set;
	      set = 0;
	    }
	  crossSections->clear();
	  delete crossSections;
	  crossSections = 0;
	}
    }

  crossSections = BuildCrossSectionsForMaterials(energyVector, energyCuts);

  if (crossSections == 0) 
    G4Exception("G4VCrossSectionHandler::BuildMeanFreePathForMaterials, crossSections = 0");

  G4VDataSetAlgorithm* algo = CreateInterpolation();
  G4VEMDataSet* materialSet = new G4CompositeEMDataSet(algo);

  G4DataVector* energies;
  G4DataVector* data;

  size_t nMaterials = G4Material::GetNumberOfMaterials();

  for (size_t m=0; m<nMaterials; m++)
    {
     energies = new G4DataVector;
      data = new G4DataVector;
      for (G4int bin=0; bin<nBins; bin++)
	{
	  G4double energy = energyVector[bin];
	  energies->push_back(energy);
	  G4VEMDataSet* matCrossSet = (*crossSections)[m];
	  G4double materialCrossSection = matCrossSet->FindValue(energy);
  
	  if (materialCrossSection > 0.)
	    {
	      data->push_back(1./materialCrossSection);
	    }
	  else
	    {
	      data->push_back(DBL_MAX);
	    }
	}
      G4VDataSetAlgorithm* algo = CreateInterpolation();
      G4VEMDataSet* dataSet = new G4EMDataSet(m,energies,data,algo,1.,1.);
      materialSet->AddComponent(dataSet);
     }

  return materialSet;
}

G4int G4VCrossSectionHandler::SelectRandomAtom(const G4Material* material, G4double e) const
{
  // Select randomly an element within the material, according to the weight 
  // determined by the cross sections in the data set
  
  G4int nElements = material->GetNumberOfElements();
  const G4ElementVector* elementVector = material->GetElementVector();

  // Special case: the material consists of one element
  if (nElements == 1) 
    {
      G4int Z = (G4int) (*elementVector)[0]->GetZ();
      return Z;
    }

  // Composite material

  G4double materialCrossSection0 = ValueForMaterial(material,e);
  //  size_t materialIndex = material->GetIndex();

  //  G4VEMDataSet* materialSet = crossSections[materialIndex];
  //  G4double materialCrossSection = materialSet->FindValue(e);

  G4double random = G4UniformRand() * materialCrossSection0;
  const G4double* nAtomsPerVolume = material->GetVecNbOfAtomsPerVolume();   
  G4double partialSumSigma = 0.;

  for ( G4int i=0 ; i < nElements ; i++ )
    { 
      G4int Z = (G4int) (*elementVector)[i]->GetZ();
      G4double crossSection = FindValue(Z,e);
      partialSumSigma += nAtomsPerVolume[i] * crossSection;
      if (random <= partialSumSigma) return Z;
    }
  // It should never get here
  return 0;
}

const G4Element* G4VCrossSectionHandler::SelectRandomElement(const G4Material* material, 
							     G4double e) const
{
  // Select randomly an element within the material, according to the weight determined
  // by the cross sections in the data set

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
      G4double materialCrossSection0 = ValueForMaterial(material,e);
      //      size_t materialIndex = material->GetIndex();
      //      G4VEMDataSet* materialSet = crossSections[materialIndex];
      //      G4double materialCrossSection = materialSet->FindValue(e);
      
      G4double random = G4UniformRand() * materialCrossSection0;
      const G4double* nAtomsPerVolume = material->GetVecNbOfAtomsPerVolume();   
      G4double partialSumSigma = 0.;
      
      for ( G4int i=0 ; i < nElements ; i++ )
	{ 
	  G4Element* element = (*elementVector)[i];
	  G4int Z = (G4int) element->GetZ();
	  G4double crossSection = FindValue(Z,e);
	  partialSumSigma += nAtomsPerVolume[i] * crossSection;
	  if (random <= partialSumSigma) return element;
	}
    }  
      // It should never end up here
      G4cout << "G4VCrossSectionHandler::SelectRandomElement - no element found" << G4endl;
      return nullElement;
}

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

  G4VEMDataSet* dataSet = 0;
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;
  pos = dataMap.find(Z);
  // The following is a workaround for STL ObjectSpace implementation, 
  // which does not support the standard and does not accept 
  // the syntax pos->first or pos->second
  // if (pos != dataMap.end()) dataSet = pos->second;
  if (pos != dataMap.end()) dataSet = (*pos).second;

  size_t nShells = dataSet->NumberOfComponents();
  for (size_t i=0; i<nShells; i++)
    {
      const G4VEMDataSet* shellDataSet = dataSet->GetComponent(i);
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

void G4VCrossSectionHandler::ActiveElements()
{
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
    G4Exception("G4VCrossSectionHandler::ActiveElements - no MaterialTable found)");
  
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

G4VDataSetAlgorithm* G4VCrossSectionHandler::CreateInterpolation()
{
  G4VDataSetAlgorithm* algorithm = new G4LogLogInterpolation;
  return algorithm;
}

G4int G4VCrossSectionHandler::NumberOfComponents(G4int Z) const
{
  G4int n = 0;
  
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;
  pos = dataMap.find(Z);
  if (pos!= dataMap.end())
    {
      G4VEMDataSet* dataSet = (*pos).second;
      n = dataSet->NumberOfComponents();
    }
  else
    {
      G4cout << "WARNING: G4VCrossSectionHandler::NumberOfComponents did not "
             << "find Z = "
             << Z << G4endl;
    }
  return n;
}



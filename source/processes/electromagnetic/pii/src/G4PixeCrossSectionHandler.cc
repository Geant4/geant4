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
// 16 Jun 2008 MGP   Created; Cross section manager for hadron impact ionization
//                   Documented in:
//                   M.G. Pia et al., PIXE Simulation With Geant4,
//                   IEEE Trans. Nucl. Sci., vol. 56, no. 6, pp. 3614-3649, Dec. 2009
//
// -------------------------------------------------------------------

#include "G4PixeCrossSectionHandler.hh"
#include "G4PhysicalConstants.hh"
#include "G4IInterpolator.hh"
#include "G4LogLogInterpolator.hh"
#include "G4IDataSet.hh"
#include "G4DataSet.hh"
#include "G4CompositeDataSet.hh"
#include "G4PixeShellDataSet.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"

#include <map>
#include <vector>
#include <fstream>
#include <sstream>


G4PixeCrossSectionHandler::G4PixeCrossSectionHandler()
{
  crossSections = 0;
  interpolation = 0;
  // Initialise with default values
  Initialise(0,"","","",1.*keV,0.1*GeV,200,MeV,barn,6,92);
  ActiveElements();
}


G4PixeCrossSectionHandler::G4PixeCrossSectionHandler(G4IInterpolator* algorithm,
						     const G4String& modelK,
						     const G4String& modelL,
						     const G4String& modelM,
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

  crossModel.push_back(modelK);
  crossModel.push_back(modelL);
  crossModel.push_back(modelM);
  
  //std::cout << "PixeCrossSectionHandler constructor - crossModel[0] = " 
  //    << crossModel[0]
  //    << std::endl;

  ActiveElements();
}

G4PixeCrossSectionHandler::~G4PixeCrossSectionHandler()
{
  delete interpolation;
  interpolation = 0;
  std::map<G4int,G4IDataSet*,std::less<G4int> >::iterator pos;

  for (pos = dataMap.begin(); pos != dataMap.end(); ++pos)
    {
      // The following is a workaround for STL ObjectSpace implementation, 
      // which does not support the standard and does not accept 
      // the syntax pos->second
      // G4IDataSet* dataSet = pos->second;
      G4IDataSet* dataSet = (*pos).second;
      delete dataSet;
    }

  if (crossSections != 0)
    {
      std::size_t n = crossSections->size();
      for (std::size_t i=0; i<n; ++i)
	{
	  delete (*crossSections)[i];
	}
      delete crossSections;
      crossSections = 0;
    }
}

void G4PixeCrossSectionHandler::Initialise(G4IInterpolator* algorithm,
					   const G4String& modelK,
					   const G4String& modelL,
					   const G4String& modelM,
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

  crossModel.push_back(modelK);
  crossModel.push_back(modelL);
  crossModel.push_back(modelM);

}

void G4PixeCrossSectionHandler::PrintData() const
{
  std::map<G4int,G4IDataSet*,std::less<G4int> >::const_iterator pos;

  for (pos = dataMap.begin(); pos != dataMap.end(); pos++)
    {
      // The following is a workaround for STL ObjectSpace implementation, 
      // which does not support the standard and does not accept 
      // the syntax pos->first or pos->second
      // G4int z = pos->first;
      // G4IDataSet* dataSet = pos->second;
      G4int z = (*pos).first;
      G4IDataSet* dataSet = (*pos).second;     
      G4cout << "---- Data set for Z = "
	     << z
	     << G4endl;
      dataSet->PrintData();
      G4cout << "--------------------------------------------------" << G4endl;
    }
}

void G4PixeCrossSectionHandler::LoadShellData(const G4String& fileName)
{
  std::size_t nZ = activeZ.size();
  for (std::size_t i=0; i<nZ; ++i)
    {
      G4int Z = (G4int) activeZ[i];     
      G4IInterpolator* algo = interpolation->Clone();
      G4IDataSet* dataSet = new G4PixeShellDataSet(Z, algo,crossModel[0],crossModel[1],crossModel[2]);

      // Degug printing
      //std::cout << "PixeCrossSectionHandler::Load - "
      //	<< Z
      //	<< ", modelK = "
      //	<< crossModel[0]
      //	<< " fileName = "
      //	<< fileName
      //	<< std::endl;

      dataSet->LoadData(fileName);
      dataMap[Z] = dataSet;
    }

  // Build cross sections for materials if not already built
  if (! crossSections)
    {
      BuildForMaterials();
    }

}

void G4PixeCrossSectionHandler::Clear()
{
  // Reset the map of data sets: remove the data sets from the map 
  std::map<G4int,G4IDataSet*,std::less<G4int> >::iterator pos;

  if(! dataMap.empty())
    {
      for (pos = dataMap.begin(); pos != dataMap.end(); ++pos)
	{
	  // The following is a workaround for STL ObjectSpace implementation, 
	  // which does not support the standard and does not accept
	  // the syntax pos->first or pos->second
	  // G4IDataSet* dataSet = pos->second;
	  G4IDataSet* dataSet = (*pos).second;
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

G4double G4PixeCrossSectionHandler::FindValue(G4int Z, G4double energy) const
{
  G4double value = 0.;
  
  std::map<G4int,G4IDataSet*,std::less<G4int> >::const_iterator pos;
  pos = dataMap.find(Z);
  if (pos!= dataMap.end())
    {
      // The following is a workaround for STL ObjectSpace implementation, 
      // which does not support the standard and does not accept 
      // the syntax pos->first or pos->second
      // G4IDataSet* dataSet = pos->second;
      G4IDataSet* dataSet = (*pos).second;
      value = dataSet->FindValue(energy);
    }
  else
    {
      G4cout << "WARNING: G4PixeCrossSectionHandler::FindValue(Z,e) did not find Z = "
	     << Z << G4endl;
    }
  return value;
}

G4double G4PixeCrossSectionHandler::FindValue(G4int Z, G4double energy, 
					      G4int shellIndex) const
{
  G4double value = 0.;

  std::map<G4int,G4IDataSet*,std::less<G4int> >::const_iterator pos;
  pos = dataMap.find(Z);
  if (pos!= dataMap.end())
    {
      // The following is a workaround for STL ObjectSpace implementation, 
      // which does not support the standard and does not accept 
      // the syntax pos->first or pos->second
      // G4IDataSet* dataSet = pos->second;
      G4IDataSet* dataSet = (*pos).second;
      if (shellIndex >= 0) 
	{
	  G4int nComponents = (G4int)dataSet->NumberOfComponents();
	  if(shellIndex < nComponents)    
	    // The value is the cross section for shell component at given energy
	    value = dataSet->GetComponent(shellIndex)->FindValue(energy);
	  else 
	    {
	      G4cout << "WARNING: G4PixeCrossSectionHandler::FindValue(Z,e,shell) did not find"
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
      G4cout << "WARNING: G4PixeCrossSectionHandler::FindValue did not find Z = "
	     << Z << G4endl;
    }
  return value;
}


G4double G4PixeCrossSectionHandler::ValueForMaterial(const G4Material* material,
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

/*
  G4IDataSet* G4PixeCrossSectionHandler::BuildMeanFreePathForMaterials(const G4DataVector*  energyCuts )
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
  std::vector<G4IDataSet*>::iterator mat;
  if (! crossSections->empty())
  {
  for (mat = crossSections->begin(); mat!= crossSections->end(); ++mat)
  {
  G4IDataSet* set = *mat;
  delete set;
  set = 0;
  }
  crossSections->clear();
  delete crossSections;
  crossSections = 0;
  }
  }

  crossSections = BuildCrossSectionsForMaterials(energyVector);

  if (crossSections == 0)
  G4Exception("G4PixeCrossSectionHandler::BuildMeanFreePathForMaterials", 
  "pii00000201",
  FatalException,
  "crossSections = 0");
  
  G4IInterpolator* algo = CreateInterpolation();
  G4IDataSet* materialSet = new G4CompositeDataSet(algo);

  G4DataVector* energies;
  G4DataVector* data;
  
  const G4ProductionCutsTable* theCoupleTable=
  G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t numOfCouples = theCoupleTable->GetTableSize();


  for (std::size_t m=0; m<numOfCouples; m++)
  {
  energies = new G4DataVector;
  data = new G4DataVector;
  for (G4int bin=0; bin<nBins; bin++)
  {
  G4double energy = energyVector[bin];
  energies->push_back(energy);
  G4IDataSet* matCrossSet = (*crossSections)[m];
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
  G4IInterpolator* algo = CreateInterpolation();
  G4IDataSet* dataSet = new G4DataSet(m,energies,data,algo,1.,1.);
  materialSet->AddComponent(dataSet);
  }

  return materialSet;
  }

*/

void G4PixeCrossSectionHandler::BuildForMaterials()
{
  // Builds a CompositeDataSet containing the mean free path for each material
  // in the material table

  G4DataVector energyVector;
  G4double dBin = std::log10(eMax/eMin) / nBins;

  for (G4int i=0; i<nBins+1; i++)
    {
      energyVector.push_back(std::pow(10., std::log10(eMin)+i*dBin));
    }

  if (crossSections != 0)
    {  // Reset the list of cross sections
      std::vector<G4IDataSet*>::iterator mat;
      if (! crossSections->empty())
	{
	  for (mat = crossSections->begin(); mat!= crossSections->end(); ++mat)
	    {
	      G4IDataSet* set = *mat;
	      delete set;
	      set = 0;
	    }
	  crossSections->clear();
	  delete crossSections;
	  crossSections = 0;
	}
    }

  crossSections = BuildCrossSectionsForMaterials(energyVector);

  if (crossSections == 0)
    G4Exception("G4PixeCrossSectionHandler::BuildForMaterials",
		"pii00000210",
		FatalException,
		", crossSections = 0");

  return;
}


G4int G4PixeCrossSectionHandler::SelectRandomAtom(const G4Material* material,
						  G4double e) const
{
  // Select randomly an element within the material, according to the weight
  // determined by the cross sections in the data set

  G4int nElements = (G4int)material->GetNumberOfElements();

  // Special case: the material consists of one element
  if (nElements == 1)
    {
      G4int Z = (G4int) material->GetZ();
      return Z;
    }

  // Composite material

  const G4ElementVector* elementVector = material->GetElementVector();
  std::size_t materialIndex = material->GetIndex();

  G4IDataSet* materialSet = (*crossSections)[materialIndex];
  G4double materialCrossSection0 = 0.0;
  G4DataVector cross;
  cross.clear();
  for ( G4int i=0; i < nElements; ++i )
    {
      G4double cr = materialSet->GetComponent(i)->FindValue(e);
      materialCrossSection0 += cr;
      cross.push_back(materialCrossSection0);
    }

  G4double random = G4UniformRand() * materialCrossSection0;

  for (G4int k=0 ; k < nElements ; ++k )
    {
      if (random <= cross[k]) return (G4int) (*elementVector)[k]->GetZ();
    }
  // It should never get here
  return 0;
}

/*
  const G4Element* G4PixeCrossSectionHandler::SelectRandomElement(const G4MaterialCutsCouple* couple,
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

  std::size_t materialIndex = couple->GetIndex();

  G4IDataSet* materialSet = (*crossSections)[materialIndex];
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
  G4cout << "G4PixeCrossSectionHandler::SelectRandomElement - no element found" << G4endl;
  return nullElement;
  }
  }
*/


G4int G4PixeCrossSectionHandler::SelectRandomShell(G4int Z, G4double e) const
{
  // Select randomly a shell, according to the weight determined by the cross sections
  // in the data set

  // Note for later improvement: it would be useful to add a cache mechanism for already
  // used shells to improve performance

  G4int shell = 0;

  G4double totCrossSection = FindValue(Z,e);
  G4double random = G4UniformRand() * totCrossSection;
  G4double partialSum = 0.;

  G4IDataSet* dataSet = 0;
  std::map<G4int,G4IDataSet*,std::less<G4int> >::const_iterator pos;
  pos = dataMap.find(Z);
  // The following is a workaround for STL ObjectSpace implementation,
  // which does not support the standard and does not accept
  // the syntax pos->first or pos->second
  // if (pos != dataMap.end()) dataSet = pos->second;
  if (pos != dataMap.end()) dataSet = (*pos).second;

  G4int nShells = (G4int)dataSet->NumberOfComponents();
  for (G4int i=0; i<nShells; ++i)
    {
      const G4IDataSet* shellDataSet = dataSet->GetComponent(i);
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

void G4PixeCrossSectionHandler::ActiveElements()
{
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
    G4Exception("G4PixeCrossSectionHandler::ActiveElements",
				  "pii00000220",
				  FatalException,
				  "no MaterialTable found");

  std::size_t nMaterials = G4Material::GetNumberOfMaterials();

  for (std::size_t mat=0; mat<nMaterials; ++mat)
    {
      const G4Material* material= (*materialTable)[mat];
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

G4IInterpolator* G4PixeCrossSectionHandler::CreateInterpolation()
{
  G4IInterpolator* algorithm = new G4LogLogInterpolator;
  return algorithm;
}

G4int G4PixeCrossSectionHandler::NumberOfComponents(G4int Z) const
{
  G4int n = 0;

  std::map<G4int,G4IDataSet*,std::less<G4int> >::const_iterator pos;
  pos = dataMap.find(Z);
  if (pos!= dataMap.end())
    {
      G4IDataSet* dataSet = (*pos).second;
      n = (G4int)dataSet->NumberOfComponents();
    }
  else
    {
      G4cout << "WARNING: G4PixeCrossSectionHandler::NumberOfComponents did not "
             << "find Z = "
             << Z << G4endl;
    }
  return n;
}


std::vector<G4IDataSet*>*
G4PixeCrossSectionHandler::BuildCrossSectionsForMaterials(const G4DataVector& energyVector)
{
  G4DataVector* energies;
  G4DataVector* data;

  std::vector<G4IDataSet*>* matCrossSections = new std::vector<G4IDataSet*>;

  //const G4ProductionCutsTable* theCoupleTable=G4ProductionCutsTable::GetProductionCutsTable();
  //std::size_t numOfCouples = theCoupleTable->GetTableSize();

  std::size_t nOfBins = energyVector.size();
  const G4IInterpolator* interpolationAlgo = CreateInterpolation();

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
    G4Exception("G4PixeCrossSectionHandler::BuildCrossSectionsForMaterials",
		"pii00000230",
		FatalException,
		"no MaterialTable found");

  std::size_t nMaterials = G4Material::GetNumberOfMaterials();

  for (std::size_t mat=0; mat<nMaterials; ++mat)
    {
      const G4Material* material = (*materialTable)[mat];
      G4int nElements = (G4int)material->GetNumberOfElements();
      const G4ElementVector* elementVector = material->GetElementVector();
      const G4double* nAtomsPerVolume = material->GetAtomicNumDensityVector();

      G4IInterpolator* algo = interpolationAlgo->Clone();

      G4IDataSet* setForMat = new G4CompositeDataSet(algo,1.,1.);

      for (G4int i=0; i<nElements; ++i) {
 
        G4int Z = (G4int) (*elementVector)[i]->GetZ();
        G4double density = nAtomsPerVolume[i];

        energies = new G4DataVector;
        data = new G4DataVector;


        for (std::size_t bin=0; bin<nOfBins; ++bin)
	  {
	    G4double e = energyVector[bin];
	    energies->push_back(e);
            G4double cross = 0.;
	    if (Z >= zMin && Z <= zMax) cross = density*FindValue(Z,e);
	    data->push_back(cross);
	  }

        G4IInterpolator* algo1 = interpolationAlgo->Clone();
        G4IDataSet* elSet = new G4DataSet(i,energies,data,algo1,1.,1.);
        setForMat->AddComponent(elSet);
      }

      matCrossSections->push_back(setForMat);
    }
  return matCrossSections;
}


G4double G4PixeCrossSectionHandler::MicroscopicCrossSection(const G4ParticleDefinition* particleDef,
							    G4double kineticEnergy,
							    G4double Z,
							    G4double deltaCut) const
{
  // Cross section formula is OK for spin=0, 1/2, 1 only !
  // Calculates the microscopic cross section in Geant4 internal units
  // Formula documented in Geant4 Phys. Ref. Manual
  // ( it is called for elements, AtomicNumber = z )

    G4double cross = 0.;

  // Particle mass and energy
    G4double particleMass = particleDef->GetPDGMass();
    G4double energy = kineticEnergy + particleMass;

  // Some kinematics
  G4double gamma = energy / particleMass;
  G4double beta2 = 1. - 1. / (gamma * gamma);
  G4double var = electron_mass_c2 / particleMass;
  G4double tMax = 2. * electron_mass_c2 * (gamma*gamma - 1.) / (1. + 2.*gamma*var + var*var);

  // Calculate the total cross section

  if ( tMax > deltaCut ) 
    {
      var = deltaCut / tMax;
      cross = (1. - var * (1. - beta2 * std::log(var))) / deltaCut;
      
      G4double spin = particleDef->GetPDGSpin() ;
      
      // +term for spin=1/2 particle
      if (spin == 0.5) 
	{
	  cross +=  0.5 * (tMax - deltaCut) / (energy*energy);
	}
      // +term for spin=1 particle
      else if (spin > 0.9 )
	{
	  cross += -std::log(var) / (3.*deltaCut) + (tMax-deltaCut) * 
	    ((5.+1./var)*0.25 /(energy*energy) - beta2 / (tMax*deltaCut))/3.;
	}
      cross *= twopi_mc2_rcl2 * Z / beta2 ;
    }

  //std::cout << "Microscopic = " << cross/barn 
  //    << ", e = " << kineticEnergy/MeV <<std:: endl; 

  return cross;
}


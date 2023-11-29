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
// 1 Aug 2001   MGP        Created
//
// 15 Jul 2009   Nicolas A. Karakatsanis
//
//                           - LoadNonLogData method was created to load only the non-logarithmic data from G4EMLOW
//                             dataset. It is essentially performing the data loading operations as in the past.
//
//                           - LoadData method was revised in order to calculate the logarithmic values of the data
//                             It retrieves the data values from the G4EMLOW data files but, then, calculates the
//                             respective log values and loads them to seperate data structures. 
//
//                           - SetLogEnergiesData method was cretaed to set logarithmic values to G4 data vectors.
//                             The EM data sets, initialized this way, contain both non-log and log values.
//                             These initialized data sets can enhance the computing performance of data interpolation
//                             operations
//
// -------------------------------------------------------------------

#include "G4CompositeEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include <fstream>
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CompositeEMDataSet::G4CompositeEMDataSet(G4VDataSetAlgorithm* argAlgorithm, 
					   G4double argUnitEnergies, 
					   G4double argUnitData, 
					   G4int argMinZ, 
					   G4int argMaxZ)
  :
  algorithm(argAlgorithm),
  unitEnergies(argUnitEnergies),
  unitData(argUnitData),
  minZ(argMinZ),
  maxZ(argMaxZ)
{
  if (algorithm == nullptr) 
  G4Exception("G4CompositeEMDataSet::G4CompositeEMDataSet",
	      "em1003",FatalException,"interpolation == 0");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CompositeEMDataSet::~G4CompositeEMDataSet()
{
  CleanUpComponents();
  if (algorithm) delete algorithm;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CompositeEMDataSet::FindValue(G4double argEnergy, G4int argComponentId) const
{
  const G4VEMDataSet* component(GetComponent(argComponentId));
 
  if (component) return component->FindValue(argEnergy);

  std::ostringstream message;
  message << "G4CompositeEMDataSet::FindValue - component " << argComponentId << " not found";
 
  G4Exception("G4CompositeEMDataSet::FindValue",
	      "em1004",FatalException,message.str().c_str());
 
  return 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CompositeEMDataSet::PrintData(void) const
{
  const G4int n = (G4int)NumberOfComponents();

  G4cout << "The data set has " << n << " components" << G4endl;
  G4cout << G4endl;
 
  G4int i(0);
 
  while (i<n)
    {
      G4cout << "--- Component " << i << " ---" << G4endl;
      GetComponent(i)->PrintData();
      ++i;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CompositeEMDataSet::SetEnergiesData(G4DataVector* argEnergies, G4DataVector* argData, G4int argComponentId)
{
  G4VEMDataSet * component(components[argComponentId]);
 
  if (component)
    {
      component->SetEnergiesData(argEnergies, argData, 0);
      return;
    }

  std::ostringstream message;
  message << "G4CompositeEMDataSet::SetEnergiesData - component " << argComponentId << " not found";
 
  G4Exception("G4CompositeEMDataSet::SetEnergiesData",
	      "em1004",FatalException,message.str().c_str());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CompositeEMDataSet::SetLogEnergiesData(G4DataVector* argEnergies, 
                                              G4DataVector* argData,
                                              G4DataVector* argLogEnergies,
                                              G4DataVector* argLogData, 
                                              G4int argComponentId)
{
  G4VEMDataSet * component(components[argComponentId]);
 
  if (component)
    {
      component->SetLogEnergiesData(argEnergies, argData, argLogEnergies, argLogData, 0);
      return;
    }

  std::ostringstream message;
  message << "G4CompositeEMDataSet::SetEnergiesData - component " << argComponentId << " not found";
 
  G4Exception("G4CompositeEMDataSet::SetLogEnergiesData",
	      "em1004",FatalException,message.str().c_str());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4CompositeEMDataSet::LoadData(const G4String& argFileName)
{
  CleanUpComponents(); 

  for (G4int z(minZ); z<maxZ; ++z)
    {
      G4VEMDataSet* component = new G4EMDataSet(z, algorithm->Clone(), unitEnergies, unitData);
      if (!component->LoadData(argFileName))
	{
	  delete component;
	  return false;
	}
      AddComponent(component);
    }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4CompositeEMDataSet::LoadNonLogData(const G4String& argFileName)
{
  CleanUpComponents(); 

  for (G4int z(minZ); z<maxZ; ++z)
    {
      G4VEMDataSet* component = new G4EMDataSet(z, algorithm->Clone(), unitEnergies, unitData);
      if (!component->LoadNonLogData(argFileName))
	{
	  delete component;
	  return false;
	}
      AddComponent(component);
    }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4CompositeEMDataSet::SaveData(const G4String& argFileName) const
{
  for (G4int z=minZ; z<maxZ; ++z)
    {
      const G4VEMDataSet* component(GetComponent(z-minZ));
  
      if (!component)
	{
	  std::ostringstream message;
	  message << "G4CompositeEMDataSet::SaveData - component " << (z-minZ) << " not found";
          G4Exception("G4CompositeEMDataSet::SaveData",
	      "em1004",FatalException,message.str().c_str());
	  return false;
	}

      if (!component->SaveData(argFileName))
	return false;
    }
 
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CompositeEMDataSet::CleanUpComponents(void)
{
  while (!components.empty())
    {
      if (components.back())
	delete components.back();
      components.pop_back();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CompositeEMDataSet::RandomSelect(G4int componentId) const
{
  G4double value = 0.;
  if (componentId >= 0 && componentId < (G4int)components.size())
    {
      const G4VEMDataSet* dataSet = GetComponent(componentId);
      value = dataSet->RandomSelect();
    }
  return value;
}

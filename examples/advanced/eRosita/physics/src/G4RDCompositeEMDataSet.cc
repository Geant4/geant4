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
// -------------------------------------------------------------------

#include "G4RDCompositeEMDataSet.hh"
#include "G4RDEMDataSet.hh"
#include "G4RDVDataSetAlgorithm.hh"
#include <fstream>
#include <sstream>

G4RDCompositeEMDataSet::G4RDCompositeEMDataSet(G4RDVDataSetAlgorithm* argAlgorithm, 
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
  if (algorithm == 0) 
    G4Exception("G4RDCompositeEMDataSet::G4RDCompositeEMDataSet()",
                "InvalidSetup", FatalException, "Interpolation == 0!");
}



G4RDCompositeEMDataSet::~G4RDCompositeEMDataSet()
{
  CleanUpComponents();
  if (algorithm) delete algorithm;
}


G4double G4RDCompositeEMDataSet::FindValue(G4double argEnergy, G4int argComponentId) const
{
  const G4RDVEMDataSet* component(GetComponent(argComponentId));
 
  if (component) return component->FindValue(argEnergy);

  std::ostringstream message;
  message << "Component " << argComponentId << " not found";
 
  G4Exception("G4RDCompositeEMDataSet::FindValue()",
              "DataNotFound", FatalException, message.str().c_str());
 
  return 0.;
}

void G4RDCompositeEMDataSet::PrintData(void) const
{
  const size_t n(NumberOfComponents());

  G4cout << "The data set has " << n << " components" << G4endl;
  G4cout << G4endl;
 
  size_t i(0);
 
  while (i<n)
    {
      G4cout << "--- Component " << i << " ---" << G4endl;
      GetComponent(i)->PrintData();
      i++;
    }
}

void G4RDCompositeEMDataSet::SetEnergiesData(G4DataVector* argEnergies, G4DataVector* argData, G4int argComponentId)
{
  G4RDVEMDataSet * component(components[argComponentId]);
 
  if (component)
    {
      component->SetEnergiesData(argEnergies, argData, 0);
      return;
    }

  std::ostringstream message;
  message << "Component " << argComponentId << " not found";
 
  G4Exception("G4RDCompositeEMDataSet::SetEnergiesData()",
              "DataNotFound", FatalException, message.str().c_str());
}

G4bool G4RDCompositeEMDataSet::LoadData(const G4String& argFileName)
{
  CleanUpComponents(); 

  for (G4int z(minZ); z<maxZ; z++)
    {
      G4RDVEMDataSet* component = new G4RDEMDataSet(z, algorithm->Clone(), unitEnergies, unitData);
      if (!component->LoadData(argFileName))
	{
	  delete component;
	  return false;
	}
      AddComponent(component);
    }
  return true;
}



G4bool G4RDCompositeEMDataSet::SaveData(const G4String& argFileName) const
{
  for (G4int z=minZ; z<maxZ; z++)
    {
      const G4RDVEMDataSet* component(GetComponent(z-minZ));
  
      if (!component)
	{
	  std::ostringstream message;
	  message << "Component " << (z-minZ) << " not found";
	  G4Exception("G4RDCompositeEMDataSet::SaveData()",
                      "DataNotFound", FatalException, message.str().c_str());
	}

      if (!component->SaveData(argFileName))
	return false;
    }
 
  return true;
}

void G4RDCompositeEMDataSet::CleanUpComponents(void)
{
  while (!components.empty())
    {
      if (components.back())
	delete components.back();
      components.pop_back();
    }
}


G4double G4RDCompositeEMDataSet::RandomSelect(G4int componentId) const
{
  G4double value = 0.;
  if (componentId >= 0 && componentId < (G4int)components.size())
    {
      const G4RDVEMDataSet* dataSet = GetComponent(componentId);
      value = dataSet->RandomSelect();
    }
  return value;
}

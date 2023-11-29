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
//  1 Aug 2001   MGP        Created
// 31 Jul 2008   MGP        Revised and renamed to G4CompositeDataSet
//
//
// -------------------------------------------------------------------

#include "G4CompositeDataSet.hh"
#include "G4DataSet.hh"
#include "G4IInterpolator.hh"
#include <fstream>
#include <sstream>

G4CompositeDataSet::G4CompositeDataSet(G4IInterpolator* algo, 
				       G4double eUnit, 
				       G4double dataUnit, 
				       G4int zMin, 
				       G4int zMax)
  :
  algorithm(algo),
  unitEnergies(eUnit),
  unitData(dataUnit),
  minZ(zMin),
  maxZ(zMax)
{
  if (algorithm == 0) 
    G4Exception("G4CompositeDataSet::G4CompositeDataSet",
		"pii00000001",
                FatalException,
		"Interpolation == 0");
}



G4CompositeDataSet::~G4CompositeDataSet()
{
  CleanUpComponents();
  if (algorithm) delete algorithm;
}


G4double G4CompositeDataSet::FindValue(G4double energy, G4int componentId) const
{
  const G4IDataSet* component(GetComponent(componentId));
 
  if (component) return component->FindValue(energy);

  std::ostringstream message;
  message << "G4CompositeDataSet::FindValue - component " << componentId << " not found";
 
   G4Exception("G4CompositeDataSet::FindValue",
	      "pii00000010",
	      FatalException,
	      message.str().c_str());
 
  return 0.;
}

void G4CompositeDataSet::PrintData(void) const
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

void G4CompositeDataSet::SetEnergiesData(G4DataVector* energies, G4DataVector* data, G4int componentId)
{
  G4IDataSet * component(components[componentId]);
 
  if (component)
    {
      component->SetEnergiesData(energies, data, 0);
      return;
    }

  std::ostringstream message;
  message << "G4CompositeDataSet::SetEnergiesData - component " << componentId << " not found";
 
  G4Exception("G4CompositeDataSet::SetEnergiesData",
	      "pii00000020",
	      FatalException,
	      message.str().c_str());

}

G4bool G4CompositeDataSet::LoadData(const G4String& argFileName)
{
  CleanUpComponents(); 

  for (G4int z(minZ); z<maxZ; ++z)
    {
      G4IDataSet* component = new G4DataSet(z, algorithm->Clone(), unitEnergies, unitData);
      if (!component->LoadData(argFileName))
	{
	  delete component;
	  return false;
	}
      AddComponent(component);
    }
  return true;
}



G4bool G4CompositeDataSet::SaveData(const G4String& argFileName) const
{
  for (G4int z=minZ; z<maxZ; ++z)
    {
      const G4IDataSet* component(GetComponent(z-minZ));
  
      if (!component)
	{
	  std::ostringstream message;
	  message << "G4CompositeDataSet::SaveData - component " << (z-minZ) << " not found";
	  G4Exception("G4CompositeDataSet::SaveData",
		      "pii00000030",
		      FatalException,
		      message.str().c_str());
	}

      if (!component->SaveData(argFileName))
	return false;
    }
 
  return true;
}

void G4CompositeDataSet::CleanUpComponents(void)
{
  while (!components.empty())
    {
      if (components.back())
	delete components.back();
      components.pop_back();
    }
}


G4double G4CompositeDataSet::RandomSelect(G4int componentId) const
{
  G4double value = 0.;
  if (componentId >= 0 && componentId < (G4int)components.size())
    {
      const G4IDataSet* dataSet = GetComponent(componentId);
      value = dataSet->RandomSelect();
    }
  return value;
}

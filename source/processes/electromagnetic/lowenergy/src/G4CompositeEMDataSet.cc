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
// $Id: G4CompositeEMDataSet.cc,v 1.9 2006/06/29 19:38:46 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 1 Aug 2001   MGP        Created
//
// -------------------------------------------------------------------

#include "G4CompositeEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include <fstream>
#include <sstream>

                                                G4CompositeEMDataSet :: G4CompositeEMDataSet(G4VDataSetAlgorithm* argAlgorithm, G4double argUnitEnergies, G4double argUnitData, G4int argMinZ, G4int argMaxZ)
:
 algorithm(argAlgorithm),
 unitEnergies(argUnitEnergies),
 unitData(argUnitData),
 minZ(argMinZ),
 maxZ(argMaxZ)
{
 if (algorithm == 0) 
  G4Exception("G4CompositeEMDataSet::G4CompositeEMDataSet - interpolation == 0");
}



                                                G4CompositeEMDataSet :: ~G4CompositeEMDataSet()
{
 CleanUpComponents();

 if (algorithm)
  delete algorithm;
}






G4double                                        G4CompositeEMDataSet :: FindValue(G4double argEnergy, G4int argComponentId) const
{
 const G4VEMDataSet * component(GetComponent(argComponentId));
 
 if (component)
  return component->FindValue(argEnergy);

 std::ostringstream message;
 message << "G4CompositeEMDataSet::FindValue - component " << argComponentId << " not found";
 
 G4Exception(message.str().c_str());
 
 return 0.;
}





void                                            G4CompositeEMDataSet :: PrintData(void) const
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





void                                            G4CompositeEMDataSet :: SetEnergiesData(G4DataVector * argEnergies, G4DataVector * argData, G4int argComponentId)
{
 G4VEMDataSet * component(components[argComponentId]);
 
 if (component)
 {
  component->SetEnergiesData(argEnergies, argData, 0);
  return;
 }

 std::ostringstream message;
 message << "G4CompositeEMDataSet::SetEnergiesData - component " << argComponentId << " not found";
 
 G4Exception(message.str().c_str());
}





G4bool                                          G4CompositeEMDataSet :: LoadData(const G4String & argFileName)
{
 CleanUpComponents(); 

 for (G4int z(minZ); z<maxZ; z++)
 {
  G4VEMDataSet * component=new G4EMDataSet(z, algorithm->Clone(), unitEnergies, unitData);
  if (!component->LoadData(argFileName))
  {
   delete component;
   return false;
  }
  
  AddComponent(component);
 }
 
 return true;
}



G4bool                                          G4CompositeEMDataSet :: SaveData(const G4String & argFileName) const
{
 for (G4int z(minZ); z<maxZ; z++)
 {
  const G4VEMDataSet * component(GetComponent(z-minZ));
  
  if (!component)
  {
   std::ostringstream message;
   message << "G4CompositeEMDataSet::SaveData - component " << (z-minZ) << " not found";
 
   G4Exception(message.str().c_str());
  }

  if (!component->SaveData(argFileName))
   return false;
 }
 
 return true;
}





void                                            G4CompositeEMDataSet :: CleanUpComponents(void)
{
 while (!components.empty())
 {
  if (components.back())
   delete components.back();

  components.pop_back();
 }
}


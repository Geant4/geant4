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
// $Id: G4TestUI.cc,v 1.1 2001-10-28 18:00:34 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 07 Oct 2001   MGP        Created
//
// -------------------------------------------------------------------

#include "G4TestUI.hh"
#include "G4ProcessTest.hh"
#include "G4Material.hh"

G4TestUI::G4TestUI(): processTest(0), polarised(false)
{
  types.push_back("compton");
  types.push_back("conversion");
  types.push_back("photoel");
  types.push_back("rayleigh");
  types.push_back("brem");
  types.push_back("ionisation");

  topics.push_back("along");
  topics.push_back("post");

  categories.push_back("lowE");
  categories.push_back("standard");
}

G4TestUI::~G4TestUI()
{
  delete processTest;
}

void G4TestUI::SelectMaterial()
{
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

 G4int nMaterials = G4Material::GetNumberOfMaterials();

  G4cout << "Select the material among the available ones: " << G4endl;
  for (G4int mat = 0; mat < nMaterials; mat++)
    {
      G4cout << mat << ") "
	     << (*theMaterialTable)[mat]->GetName()
	     << G4endl;
    }
  G4int materialId;
  G4cin >> materialId;

  G4Material* material = (*theMaterialTable)[materialId] ;

  G4cout << "The selected material is: " << material->GetName() << G4endl;
}

void G4TestUI::SelectProcess()
{
  selectType();
  selectCategory();
  isPolarised();
}

void G4TestUI::SelectProcessType()
{
  G4cout << "Process to be tested: " << G4endl
	 << "Compton [1], GammaConversion [2], Photoelectric [3], Rayleigh [4]" 
	 << G4endl
	 << "Bremsstrahlung [5], eIonisation [6]" << G4endl;
  G4cin >> processType;
  if (processType < 1 || processType > 6) G4Exception("Wrong input");

  type = processType - 1;
}

void G4TestUI::SelectProcessCategory()
{
  G4int selection;
  G4cout << "LowEnergy [1] or Standard [2]" << G4endl;
  G4cin >> processSelection;
  if (processSelection < 1 || processSelection > 2) G4Exception("Wrong input");

  category = selection - 1;
}

void G4TestUI::IsPolarised()
{
  if (type < 5)
    {
      G4int isPolarised;
      G4cout << "Polarised processes are available for: Compton" 
	     << "Not Polarised [0] or Polarised [1] Process?"	       
	     << G4endl;
      G4cin >> isPolarised;
      if (isPolarised > 0) polarised = true;
    }
}

void G4TestUI::SelectTestTopic()
{
  G4int iStep;
  G4cout << "PostStep [1] or AlongStep [2] test?" << G4endl;
  G4cin >> iStep;
  topic = iStep - 1;
}

void G4TestUI::SelectNumberOfIterations()
{
  G4cout << "How many iterations? " << G4endl;
  G4cin >> nIterations;
  if (nIterations <= 0) G4Exception("Wrong input");
}

const G4Material* G4TestUI::GetSelectedMaterial()
{
  G4Material* material = (*theMaterialTable)[materialId] ;
  return material;
}

G4ProcessTest* G4TestUI::GetSelectedTest()
{
  G4ProcessTest* test = new G4ProcessTest;
  G4String  = categories[category];
  
  if (type == 0)
    {
      test = new G4ComptonTest(cat,isPolarised);
    }
  else if (type == 1)
    {
      test = new G4GammaConversionTest(cat,isPolarised);
    }
  else if (type == 2)
    {
      test = new G4PhotoelectricTest(cat,isPolarised);
    }
  else if (type == 3)
    {
      test = new G4RayleighTest(cat,isPolarised);
    }
  else if (type == 4)
    {
      test = new G4BremsstralungTest(cat);
    }
  else if (type == 5)
    {
      test = new G4eIonisationTest(cat);
    }

  return test;
}

G4int G4TestUI::GetNIterations()
{
  return nIterations;
}

const G4String& G4TestUI::GetTestTopic()
{
  return topics[topic];
}

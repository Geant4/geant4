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
// $Id: G4TestUI.cc,v 1.4 2001-10-30 08:35:52 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 07 Oct 2001   MGP        Created
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4TestUI.hh"
#include "G4ProcessTest.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"

G4TestUI::G4TestUI(): nIterations(0), polarised(false)
{
  types.push_back("compton");
  types.push_back("conversion");
  types.push_back("photoelectric");
  types.push_back("rayleigh");
  types.push_back("bremsstrahlung");
  types.push_back("ionisation");

  topics.push_back("along");
  topics.push_back("post");

  categories.push_back("lowE");
  categories.push_back("standard");
}

G4TestUI::~G4TestUI()
{ }

void G4TestUI::configure()
{
  selectNumberOfIterations();
  selectProcess();
  selectMaterial();
  selectEnergyRange();
  selectTestTopic(); 
}

void G4TestUI::selectMaterial()
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
  G4cin >> materialId;

  G4Material* material = (*theMaterialTable)[materialId] ;
  G4cout << "The selected material is: " << material->GetName() << G4endl;
}

void G4TestUI::selectProcess()
{
  selectProcessType();
  selectProcessCategory();
  isPolarised();
}

void G4TestUI::selectProcessType()
{
  G4int processType;
  G4cout << "Process to be tested: " << G4endl
	 << "Compton [1], GammaConversion [2], Photoelectric [3], Rayleigh [4]" 
	 << G4endl
	 << "Bremsstrahlung [5], eIonisation [6]" << G4endl;
  G4cin >> processType;
  if (processType < 1 || processType > 6) G4Exception("Wrong input");

  type = processType - 1;
}

void G4TestUI::selectProcessCategory()
{
  G4int selection;
  G4cout << "LowEnergy [1] or Standard [2]" << G4endl;
  G4cin >> selection;
  if (selection < 1 || selection > 2) G4Exception("Wrong input");

  category = selection - 1;
}

void G4TestUI::isPolarised()
{
  if (type == 0)
    {
      G4int isPolarised;
      G4cout << "Polarised processes are available for: Compton" 
	     << G4endl
	     << "Not Polarised [0] or Polarised [1] Process?"	       
	     << G4endl;
      G4cin >> isPolarised;
      if (isPolarised > 0) polarised = true;
    }
}

void G4TestUI::selectTestTopic()
{
  G4int iStep;
  G4cout << "PostStep [1] or AlongStep [2] test?" << G4endl;
  G4cin >> iStep;
  topic = iStep - 1;
}

void G4TestUI::selectNumberOfIterations()
{
  G4cout << "How many iterations? " << G4endl;
  G4cin >> nIterations;
  if (nIterations <= 0) G4Exception("Wrong input");
}

void G4TestUI::selectEnergyRange()
{
  G4cout << "Select min and max energy (MeV)" << G4endl;
  G4cin >> eMin >> eMax;
  if (eMin <= 0. || eMax < eMin) G4Exception("Wrong input");
  eMin = eMin * MeV;
  eMax = eMax * MeV;
}

const G4Material* G4TestUI::getSelectedMaterial() const
{
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  const G4Material* material = (*theMaterialTable)[materialId] ;

  return material;
}

G4int G4TestUI::getNumberOfIterations() const
{
  return nIterations;
}

const G4String& G4TestUI::getTestTopic() const 
{
  return topics[topic];
}

const G4String& G4TestUI::getProcessType() const 
{
  return types[type];
}

const G4String& G4TestUI::getProcessCategory() const
{
  return categories[category];
}

G4bool G4TestUI::getPolarisationSelection() const
{
  return polarised;
}

G4ParticleDefinition* G4TestUI::getParticleDefinition() const 
{
  G4ParticleDefinition* def = 0;
  if (type <= 4)
    {
      def = G4Gamma::GammaDefinition();
    }
  else
    {
      def = G4Electron::ElectronDefinition();
    }
  return def;
}



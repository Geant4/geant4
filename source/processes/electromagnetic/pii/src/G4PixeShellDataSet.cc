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
// 09.10.01   V.Ivanchenko Add case z=0
//  9 Mar 2008   MGP        Cleaned up unreadable code modified by former developer
//                          (Further clean-up needed) 
// 31 Jul 2008   MGP        Revised 
//
// -------------------------------------------------------------------

#include "G4PixeShellDataSet.hh"
#include "G4DataSet.hh"
#include "G4IInterpolator.hh"
#include <fstream>
#include <sstream>


G4PixeShellDataSet::G4PixeShellDataSet(G4int zeta, 
				       G4IInterpolator* algo, 
				       const G4String& modelK,
				       const G4String& modelL,
				       const G4String& modelM,
				       G4double eUnit, 
				       G4double dataUnit):
  z(zeta),
  algorithm(algo),
  unitEnergies(eUnit),
  unitData(dataUnit)
{
  if (algorithm == 0) G4Exception("G4PixeShellDataSet::G4PixeShellDataSet",
				  "pii00000301",
				  FatalException,
				  "interpolation == 0");

  crossModel.push_back(modelK);
  crossModel.push_back(modelL);
  crossModel.push_back(modelM);

  shellName.push_back("k");
  shellName.push_back("l");
  shellName.push_back("m");

  std::size_t sizeK = modelK.size();
  std::size_t sizeL = modelL.size();
  std::size_t sizeM = modelM.size();
  
  if (sizeK > 0) subShellName.push_back("k");

  if (sizeK > 0 && sizeL > 0)
    {
      subShellName.push_back("l1");
      subShellName.push_back("l2");
      subShellName.push_back("l3");
    }
  if (sizeK > 0 && sizeL > 0 && sizeM >0)
    {
      subShellName.push_back("m1");
      subShellName.push_back("m2");
      subShellName.push_back("m3");
      subShellName.push_back("m4");
      subShellName.push_back("m5");
    }
}


G4PixeShellDataSet::~G4PixeShellDataSet()
{
  CleanUpComponents();
  if (algorithm) delete algorithm;
}


G4double G4PixeShellDataSet::FindValue(G4double energy, G4int /* componentId */) const
{
  // Returns the sum over the shells corresponding to e
  G4double value = 0.;

  std::vector<G4IDataSet *>::const_iterator i(components.begin());
  std::vector<G4IDataSet *>::const_iterator end(components.end());

  while (i != end)
    {
      value += (*i)->FindValue(energy);
      i++;
    }
  return value;
}


void G4PixeShellDataSet::PrintData(void) const
{
  const G4int n = (G4int)NumberOfComponents();

  G4cout << "The data set has " << n << " components" << G4endl;
  G4cout << G4endl;
 
  G4int i = 0;
 
  while (i < n)
    {
      G4cout << "--- Component " << i << " ---" << G4endl;
      GetComponent(i)->PrintData();
      ++i;
    }
}


void G4PixeShellDataSet::SetEnergiesData(G4DataVector* energies, 
					 G4DataVector* data, 
					 G4int componentId)
{
  G4IDataSet* component = components[componentId];
 
  if (component)
    {
      component->SetEnergiesData(energies, data, 0);
      return;
    }

  std::ostringstream message;
  message << "G4PixeShellDataSet::SetEnergiesData - component " << componentId << " not found";
 
  G4Exception("G4PixeShellDataSet::SetEnergiesData",
	      "pii000000310",
	      FatalException,
	      message.str().c_str());
}


G4bool G4PixeShellDataSet::LoadData(const G4String& file)
{
  CleanUpComponents();

  // Load shell cross sections
  
  std::size_t nShells = subShellName.size();
  
  for (std::size_t subShellIndex=0; subShellIndex<nShells; ++subShellIndex)
    {
      G4String subName = subShellName[subShellIndex];    
      G4String fullFileName = FullFileName(file,subName);

      // Create component DataSet with the data from the current subshell
      G4IDataSet* dataSet = new G4DataSet(z,algorithm);
      dataSet->LoadData(fullFileName);
      
      // Add component to the ShellDataSet
      AddComponent(dataSet);
    }

  return true;
}


G4bool G4PixeShellDataSet::SaveData(const G4String& /* file */) const
{
  // Dummy implementation
  return true;
}


void G4PixeShellDataSet::CleanUpComponents(void)
{
  while (!components.empty())
    {
      if (components.back()) delete components.back();
      components.pop_back();
    }
}


G4String G4PixeShellDataSet::FullFileName(const G4String& file,
					  const G4String& subShell) const
{
  const char* path = G4FindDataDir("G4PIIDATA");
  if (!path)
    G4Exception("G4PixeShellDataSet::FullFileName",
				  "pii00000320",
				  FatalException,
				  "G4PIIDATA environment variable not set");
  
  // Identify the shell this subshell belongs to
  G4int shellIndex = TranslateShell(subShell);
  G4String shellString = shellName[shellIndex];
  G4String shellModel = crossModel[shellIndex];

  std::ostringstream fullFileName;
 
  fullFileName 
    //<< path 
	       << "pixe/" 
	       << file
	       << '/' 
	       << shellString
	       << '/'
	       << shellModel
	       << '/'
	       << subShell
	       << '-' ;
//	       << z 
	//       << ".dat";
                    
  G4String test(fullFileName.str().c_str());
  // std::cout << "PixeShellDataSet - Reading data from file " << test << std::endl;

  return G4String(fullFileName.str().c_str());
}

G4int G4PixeShellDataSet::TranslateShell(const G4String& subShell) const
{
  // By default return K shell
  G4int index = 0;

  if (subShell == "l1" || subShell == "l2" || subShell == "l3" ) index = 1;
  if (subShell == "m1" || 
      subShell == "m2" ||
      subShell == "m3" ||
      subShell == "m4" || 
      subShell == "m5" ) index = 2;
  return index;
}

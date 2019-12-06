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
//
// -------------------------------------------------------------------

#include "G4RDShellEMDataSet.hh"
#include "G4RDEMDataSet.hh"
#include "G4RDVDataSetAlgorithm.hh"
#include <fstream>
#include <sstream>


G4RDShellEMDataSet::G4RDShellEMDataSet(G4int zeta, G4RDVDataSetAlgorithm* algo, 
				   G4double eUnit, 
				   G4double dataUnit)
  :
  z(zeta),
  algorithm(algo),
  unitEnergies(eUnit),
  unitData(dataUnit)
{
  if (algorithm == 0)
    G4Exception("G4RDShellEMDataSet::G4RDShellEMDataSet()", "InvalidSetup",
                FatalException, "Interpolation == 0!");
}


G4RDShellEMDataSet::~G4RDShellEMDataSet()
{
  CleanUpComponents();
  if (algorithm) delete algorithm;
}


G4double G4RDShellEMDataSet::FindValue(G4double energy, G4int /* componentId */) const
{
  // Returns the sum over the shells corresponding to e
  G4double value = 0.;

  std::vector<G4RDVEMDataSet *>::const_iterator i(components.begin());
  std::vector<G4RDVEMDataSet *>::const_iterator end(components.end());

  while (i != end)
    {
      value += (*i)->FindValue(energy);
      i++;
    }

  return value;
}


void G4RDShellEMDataSet::PrintData(void) const
{
  const size_t n = NumberOfComponents();

  G4cout << "The data set has " << n << " components" << G4endl;
  G4cout << G4endl;
 
  size_t i = 0;
 
  while (i < n)
    {
      G4cout << "--- Component " << i << " ---" << G4endl;
      GetComponent(i)->PrintData();
      i++;
    }
}


void G4RDShellEMDataSet::SetEnergiesData(G4DataVector* energies, 
				       G4DataVector* data, 
				       G4int componentId)
{
  G4RDVEMDataSet* component = components[componentId];
 
  if (component)
    {
      component->SetEnergiesData(energies, data, 0);
      return;
    }

  std::ostringstream message;
  message << "Component " << componentId << " not found";
 
  G4Exception("G4RDShellEMDataSet::SetEnergiesData()", "DataNotFound",
              FatalException, message.str().c_str());
}


G4bool G4RDShellEMDataSet::LoadData(const G4String& file)
{
  CleanUpComponents();

  G4String fullFileName = FullFileName(file);
  std::ifstream in(fullFileName);

  if (!in.is_open())
    {
      G4String message("Data file \"");
      message += fullFileName;
      message += "\" not found";
      G4Exception("G4RDShellEMDataSet::LoadData()", "DataNotFound",
                  FatalException, message);
    }

  G4DataVector* energies = 0;
  G4DataVector* data = 0;

  G4double a = 0.;
  G4int shellIndex = 0;
  bool energyColumn = true;

  do
    {
      in >> a;
  
      if (a == -1)
	{
	  if (energyColumn && energies!=0)
	    {
	      AddComponent(new G4RDEMDataSet(shellIndex, energies, data, algorithm->Clone(), unitEnergies, unitData));
	      energies = 0;
	      data = 0;
	    }
   
	  energyColumn = (!energyColumn);
	}
      else if (a != -2)
	{
	  if (energies == 0)
	    {
	      energies = new G4DataVector;
	      data = new G4DataVector;
	    }
  
	  if (energyColumn)
	    energies->push_back(a * unitEnergies);
	  else
	    data->push_back(a * unitData);

	  energyColumn = (!energyColumn);
	}
    }
  while (a != -2);

  return true;
}


G4bool G4RDShellEMDataSet::SaveData(const G4String& file) const
{
  G4String fullFileName = FullFileName(file);
  std::ofstream out(fullFileName);

  if (!out.is_open())
    {
      G4String message("Cannot open \"");
      message += fullFileName;
      message += "\"";
      G4Exception("G4RDEMDataSet::SaveData()", "CannotOpenFile",
                  FatalException, message);
    }
 
  const size_t n = NumberOfComponents();
  size_t k = 0;
 
  while (k < n)
    {
      const G4RDVEMDataSet* component = GetComponent(k);
  
      if (component)
	{
	  const G4DataVector& energies = component->GetEnergies(0);
	  const G4DataVector& data = component->GetData(0);
 
	  G4DataVector::const_iterator i = energies.begin();
	  G4DataVector::const_iterator endI = energies.end();
	  G4DataVector::const_iterator j = data.begin();
  
	  while (i != endI)
	    {
	      out.precision(10);
	      out.width(15);
	      out.setf(std::ofstream::left);
	      out << ((*i)/unitEnergies) << ' ';

	      out.precision(10);
	      out.width(15);
	      out.setf(std::ofstream::left);
	      out << ((*j)/unitData) << std::endl;
	      i++;
	      j++;
	    }
	}
  
      out.precision(10);
      out.width(15);
      out.setf(std::ofstream::left);
      out << -1.f << ' ';

      out.precision(10);
      out.width(15);
      out.setf(std::ofstream::left);
      out << -1.f << std::endl;
  
      k++;
    }
 
  out.precision(10);
  out.width(15);
  out.setf(std::ofstream::left);
  out << -2.f << ' ';

  out.precision(10);
  out.width(15);
  out.setf(std::ofstream::left);
  out << -2.f << std::endl;

  return true;
}


void G4RDShellEMDataSet::CleanUpComponents(void)
{
  while (!components.empty())
    {
      if (components.back()) delete components.back();
      components.pop_back();
    }
}


G4String G4RDShellEMDataSet::FullFileName(const G4String& fileName) const
{
  char* path = std::getenv("G4LEDATA");
  if (!path)
    G4Exception("G4RDShellEMDataSet::FullFileName()", "InvalidSetup",
                FatalException, "G4LEDATA environment variable not set!");
  
  std::ostringstream fullFileName;
 
  fullFileName << path << '/' << fileName << z << ".dat";
                      
  return G4String(fullFileName.str().c_str());
}

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
// $Id: G4ShellEMDataSet.cc,v 1.8 2002-05-28 09:20:21 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 1 Aug 2001   MGP        Created
// 09.10.01   V.Ivanchenko Add case z=0
//
// -------------------------------------------------------------------

#include "G4ShellEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "g4std/fstream"
#include "g4std/strstream"


G4ShellEMDataSet::G4ShellEMDataSet(G4int Z,
				   const G4VDataSetAlgorithm* interpolation,
				   G4double unitE, G4double unitData)
  :z(Z), algorithm(interpolation)
{
  nComponents = 0;
  unit1 = unitE;
  unit2 = unitData;
}

G4ShellEMDataSet::G4ShellEMDataSet(G4int Z, const G4String& dataFile,
				   const G4VDataSetAlgorithm* interpolation,
				   G4double unitE, G4double unitData)
  :z(Z), algorithm(interpolation)
{
  nComponents = 0;
  unit1 = unitE;
  unit2 = unitData;
  LoadData(dataFile);
}

G4ShellEMDataSet::~G4ShellEMDataSet()
{ 
  for (size_t i=0; i<nComponents; i++)
    {
      G4VEMDataSet* dataSet = components[i];
      delete dataSet;
    }
  delete algorithm;
}

G4double G4ShellEMDataSet::FindValue(G4double e, G4int id) const
{
  // Returns the sum over the shells corresponding to e
  G4double value = 0.;

  for (size_t i=0; i<nComponents; i++)
    {
      G4VEMDataSet* component = components[i];
      G4double shellValue = component->FindValue(e);
      value = value + shellValue;
    }

  return value;
}

void G4ShellEMDataSet::PrintData() const
{
  G4cout << "The data set has " << nComponents << " components" << G4endl;

  for (size_t i=0; i<nComponents; i++)
  {
    G4cout << "--- Component " << i << " ---" << G4endl;
    G4VEMDataSet* component = components[i];
    component->PrintData();
  }
}

void G4ShellEMDataSet::LoadData(const G4String& fileName)
{ 
  // Build the complete string identifying the file with the data set
  
  char nameChar[100] = {""};
  G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
  
  if (z != 0)  ost << fileName << z << ".dat";
  else   ost << fileName << ".dat";

  G4String name(nameChar);
  
  char* path = getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep("G4ShellEMDataSet - G4LEDATA environment variable not set");
      G4Exception(excep);
    }
  
  G4String pathString(path);
  G4String separator("/" );
  G4String dirFile = pathString + separator + name;
  G4std::ifstream file(dirFile);
  G4std::filebuf* lsdp = file.rdbuf();

  if (! (lsdp->is_open()) )
    {
      G4String s1("G4ShellEMDataSet - data file: ");
      G4String s2(" not found");
      G4String excep = s1 + dirFile + s2;
      G4Exception(excep);
    }

  G4double a = 0;
  G4int k = 1;
  G4int s = 0;
  
  G4int shellIndex = 0;
  G4DataVector* energies = new G4DataVector;
  G4DataVector* data = new G4DataVector;

  do {
    file >> a;
    G4int nColumns = 2;
    if (a == -1)
      {
	if (s == 0)
	  {
	    // End of a shell data set
	    G4VDataSetAlgorithm* algo = algorithm->Clone();
	    G4VEMDataSet* dataSet = new G4EMDataSet(shellIndex,energies,data,algo);
	    AddComponent(dataSet);
	    // Start of new shell data set
	    energies = new G4DataVector;
            data = new G4DataVector;
            shellIndex++;	    
	  }      
	s++;
	if (s == nColumns)
	{
	  s = 0;
	}
      }
    else if (a == -2)
      {
	// End of file; delete the empty vectors created when encountering the last -1 -1 row
	delete energies;
	delete data;
      }
    else
      {
	// 1st column is energy
	if(k%nColumns != 0)
	  {	    
	    G4double e = a * unit1;
	    energies->push_back(e);
	    k++;
	  }
	else if (k%nColumns == 0)
	  {
	    // 2nd column is cross section 
	    G4double value = a * unit2;
	    data->push_back(value);
	    k = 1;
	  }
      }
  } while (a != -2); // end of file
  file.close();
}

void G4ShellEMDataSet::AddComponent(G4VEMDataSet* component)
{ 
  components.push_back(component);
  nComponents++;
}

const G4DataVector& G4ShellEMDataSet::GetEnergies(G4int i) const
{
  const G4VEMDataSet* component = GetComponent(i);
  return (component->GetEnergies(i));
}

const G4DataVector& G4ShellEMDataSet::GetData(G4int i) const 
{
  const G4VEMDataSet* component = GetComponent(i);
  return (component->GetData(i));
}

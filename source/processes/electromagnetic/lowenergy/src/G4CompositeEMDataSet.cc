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
// $Id: G4CompositeEMDataSet.cc,v 1.6 2002-05-28 09:20:18 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "g4std/fstream"
#include "g4std/strstream"


G4CompositeEMDataSet::G4CompositeEMDataSet(G4VDataSetAlgorithm* interpolation,
					   G4double unitE, G4double unitData,
					   G4int minZ, G4int maxZ)
  :algorithm(interpolation), unit1(unitE), unit2(unitData), zMin(minZ), zMax(maxZ)
{
  nComponents = 0;
}

G4CompositeEMDataSet::G4CompositeEMDataSet(const G4String& dataFile,
					   G4VDataSetAlgorithm* interpolation,
					   G4double unitE, G4double unitData,
					   G4int minZ, G4int maxZ)
  : algorithm(interpolation), unit1(unitE), unit2(unitData), zMin(minZ), zMax(maxZ)
{
  nComponents = 0;
  LoadData(dataFile);
}

G4CompositeEMDataSet::~G4CompositeEMDataSet()
{ 
  for (size_t i=0; i<nComponents; i++)
    {
      G4VEMDataSet* dataSet = components[i];
      delete dataSet;
    }
  delete algorithm;
}

G4double G4CompositeEMDataSet::FindValue(G4double e, G4int id) const
{
  // Returns the value in component id corresponding to e
  G4double value = 0.;

  G4VEMDataSet* component = components[id];
  if (component != 0)
    {
      value = component->FindValue(e);
    }
  else
    {
      G4cout << "WARNING - G4CompositeEMDataSet::FindValue - component "
	     << id << " not found" << G4endl; 
    }

  return value;
}

void G4CompositeEMDataSet::PrintData() const
{
  G4cout << "The data set has " << nComponents << " components" << G4endl;

  for (size_t i=0; i<nComponents; i++)
  {
    G4cout << "--- Component " << i << " ---" << G4endl;
    G4VEMDataSet* component = components[i];
    component->PrintData();
  }
}

void G4CompositeEMDataSet::LoadData(const G4String& fileName)
{ 

  for (G4int Z=zMin; Z<zMax; Z++)
    {
      // Build the complete string identifying the file with the data set
      
      char nameChar[100] = {""};
      G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
      
      ost << fileName << Z << ".dat";
      
      G4String name(nameChar);
      
      char* path = getenv("G4LEDATA");
      if (!path)
	{ 
	  G4String excep("G4CompositeEMDataSet - G4LEDATA environment variable not set");
	  G4Exception(excep);
	}
      
      G4String pathString(path);
      G4String separator = G4String("/");
      G4String dirFile = pathString + separator + name;
      G4std::ifstream file(dirFile);
      G4std::filebuf* lsdp = file.rdbuf();
      
      if (! (lsdp->is_open()) )
	{
	  G4String s1("G4CompositeEMDataSet - data file: ");
          G4String s2(" not found");
	  G4String excep = s1 + dirFile + s2;
	  G4Exception(excep);
	}
      G4double a = 0;
      G4int k = 1;
      G4DataVector* energies = new G4DataVector;
      G4DataVector* data = new G4DataVector;
      do
	{
	  file >> a;
	  G4int nColumns = 2;
	  // The file is organized into two columns:
	  // 1st column is the energy
	  // 2nd column is the corresponding value
	  // The file terminates with the pattern: -1   -1
	  //                                       -2   -2
	  if (a == -1 || a == -2)
	    {
	    }
	  else
	    {
	      if (k%nColumns != 0)
		{	
		  G4double e = a * unit1;
		  energies->push_back(e);
		  k++;
		}
	      else if (k%nColumns == 0)
		{
		  G4double value = a * unit2;
		  data->push_back(value);
		  k = 1;
		}
	    }
	} while (a != -2); // end of file
      
      file.close();

      G4VDataSetAlgorithm* algo = algorithm->Clone();
      G4VEMDataSet* dataSet = new G4EMDataSet(Z,energies,data,algo);
      AddComponent(dataSet);
    }
}

void G4CompositeEMDataSet::AddComponent(G4VEMDataSet* component)
{ 
  components.push_back(component);
  nComponents++;
}

const G4DataVector& G4CompositeEMDataSet::GetEnergies(G4int i) const
{
  const G4VEMDataSet* component = GetComponent(i);
  return (component->GetEnergies(i));
}

const G4DataVector& G4CompositeEMDataSet::GetData(G4int i) const 
{
  const G4VEMDataSet* component = GetComponent(i);
  return (component->GetData(i));
}

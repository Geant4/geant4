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
// $Id: G4BremsstrahlungParameters.cc,v 1.1 2001-08-20 16:37:37 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
//
// -------------------------------------------------------------------

#include "G4BremsstrahlungParameters.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "g4std/fstream"
#include "g4std/strstream"

G4BremsstrahlungParameters:: G4BremsstrahlungParameters(G4int minZ, G4int maxZ)
  : zMin(minZ), zMax(maxZ)
{
  LoadData();
}

G4BremsstrahlungParameters::~G4BremsstrahlungParameters()
{ }

G4double G4BremsstrahlungParameters::ParameterA(G4int Z, G4double energy) const
{
  // To be implemented
  G4double value = 0.;
  return value;
}

G4double G4BremsstrahlungParameters::ParameterB(G4int Z, G4int id) const
{
  if (Z < (zMin-1) || Z > (zMax-1))
    G4Exception("G4BremsstrahlungParameters::FindValue - index outside boundaries");

  G4double value;

  if (id == 1) 
    {
      value = paramB0[Z-1];
    }
  else
    {
      value = paramB1[Z-1];
    }
  return value;
}



void G4BremsstrahlungParameters::LoadData()
{
  // Partial implementation, to be completed

  // Build the complete string identifying the file with the data set
  
  char nameChar[100] = {""};
  G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
  
  G4String fileName = "brem/br-co-b";
  ost << fileName  << ".dat";
  
  G4String name(nameChar);
  
  char* path = getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep = "G4BremsstrahlungParameters - G4LEDATA environment variable not set";
      G4Exception(excep);
    }
  
  G4String pathString(path);
  G4String dirFile = pathString + "/" + name;
  G4std::ifstream file(dirFile);
  G4std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
	{
	  G4String excep = "G4BremsstrahlungParameters - data file: " + dirFile + " not found";
	  G4Exception(excep);
	}

  G4int nParam = 2;
  G4double a = 0;
  G4int k = 1;

  do
    {
      file >> a;
      // The file is organized into two columns:
      // 1st column is the energy
      // 2nd column is the corresponding value
      // The file terminates with the pattern: -1   -1
      //                                       -2   -2

      if (k%nParam != 0)
	{	
	  paramB0.push_back(a);
	  k++;
	}
      else if (k%nParam == 0)
	{
	  paramB1.push_back(a);
	  k = 1;
	}
    } while (a != -2); // end of file
  
  file.close();
}

void G4BremsstrahlungParameters::PrintData() const
{
  // To be completed
  for (G4int i=zMin-1; i<zMax-1; i++)
    {
      G4double b0 = paramB0[i];
      G4double b1 = paramB1[i] ;
      G4cout << "Z = " << i << " - Parameter B: "
	     << b0
	     << ",  "
	     << b1
	     << G4endl; 
    }
}

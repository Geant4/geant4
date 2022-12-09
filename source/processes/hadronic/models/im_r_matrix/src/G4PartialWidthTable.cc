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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//      For information related to this code contact:
//
//      File name:     G4PartialWidthTable
//
//      Author:        
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicException.hh"
#include "G4ios.hh"
#include "G4PartialWidthTable.hh"
#include "G4KineticTrack.hh"

G4PartialWidthTable::G4PartialWidthTable(const G4double* energies, G4int entries) : nEnergies(entries)
{ 
  // Fill the vector with tabulated energies
  G4int i;
  for (i=0; i<entries; i++)
    {
      G4double e = *(energies + i) * GeV;
      energy.push_back(e);
    }
}


G4PartialWidthTable::~G4PartialWidthTable()
{ }


G4bool G4PartialWidthTable::operator==(const G4PartialWidthTable &right) const
{
  return (this == (G4PartialWidthTable *) &right);
}


G4bool G4PartialWidthTable::operator!=(const G4PartialWidthTable &right) const
{
  return (this != (G4PartialWidthTable *) &right);
}


G4int G4PartialWidthTable::NumberOfChannels() const
{
  return (G4int)widths.size();
}


const G4PhysicsVector* G4PartialWidthTable::Width(const G4String& name1, const G4String& name2) const
{
  // Returned pointer is not owned by the client
  G4PhysicsVector* width = 0;
  std::size_t n = 0;
  std::size_t entries = widths.size();
  for (std::size_t i=0; i<entries; ++i)
    {
      if ( (daughter1[i] == name1 && daughter2[i] == name2) ||
	   (daughter2[i] == name1 && daughter1[i] == name2) )
	{
	  width = (G4PhysicsVector*) (widths[i]);
	  ++n;
	}
    }
  if (n > 1) throw G4HadronicException(__FILE__, __LINE__, "G4PartialWidthTable::Width - ambiguity");

  return width;
}


void G4PartialWidthTable::AddWidths(const G4double* channelWidth, 
				    const G4String& name1, const G4String& name2)
{
  G4PhysicsFreeVector* width = new G4PhysicsFreeVector(nEnergies);
  G4int i;
  for (i=0; i<nEnergies; i++)
    {
      G4double value = *(channelWidth+ i) * GeV;
      G4double e = energy[i];
      width->PutValue(i,e,value);
    }	            

  widths.push_back(width);
  daughter1.push_back(name1);
  daughter2.push_back(name2);

  return;
}


void G4PartialWidthTable::Dump() const
{
  std::size_t entries = widths.size();

  for (std::size_t i=0; i<entries; ++i)
    {
      G4cout << " Channel " << i  << ": " 
             << daughter1[i] << " " << daughter2[i] << G4endl;
      G4PhysicsFreeVector* width = widths[i];
      for (G4int j=0; j<nEnergies; ++j)
	{
	  G4bool dummy = false;
	  G4double e = energy[i];
	  G4double w = width->GetValue(e,dummy);
	  G4cout << j << ") Energy = " << e << ", Width = " << w << G4endl;
	}
    }
  return;
}

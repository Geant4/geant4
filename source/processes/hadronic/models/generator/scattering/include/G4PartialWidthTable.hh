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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:     G4PartialWidthTable
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
//  
//
// -------------------------------------------------------------------

#ifndef G4PARTIALWIDTHTABLE_HH
#define G4PARTIALWIDTHTABLE_HH

#include "globals.hh"
#include "G4PhysicsFreeVector.hh"
#include "g4std/vector"

class G4PartialWidthTable 
{

public:

  G4PartialWidthTable(const G4double* energies, G4int entries);

  virtual ~G4PartialWidthTable();

  G4bool operator==(const G4PartialWidthTable &right) const;
  G4bool operator!=(const G4PartialWidthTable &right) const;

  G4int NumberOfChannels() const;
 
  // Ownership of the pointer is not transferred to the client
  const G4PhysicsVector* Width(const G4String& name1, const G4String& name2) const;

  void AddWidths(const G4double* widths, const G4String& name1, const G4String& name2);

  void Dump() const;


protected:

private:  

  G4PartialWidthTable(const G4PartialWidthTable &right);
  const G4PartialWidthTable& operator=(const G4PartialWidthTable &right);

  G4int nEnergies;
  G4std::vector<G4double> energy;
  G4std::vector<G4PhysicsFreeVector*> widths;
  G4std::vector<G4String> daughter1;
  G4std::vector<G4String> daughter2;

};

#endif



















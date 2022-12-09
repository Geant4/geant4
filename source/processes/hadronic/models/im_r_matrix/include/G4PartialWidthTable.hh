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
#include <vector>

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
  std::vector<G4double> energy;
  std::vector<G4PhysicsFreeVector*> widths;
  std::vector<G4String> daughter1;
  std::vector<G4String> daughter2;

};

#endif

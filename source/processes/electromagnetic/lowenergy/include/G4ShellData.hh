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
// $Id: G4ShellData.hh,v 1.1 2001-08-20 16:36:01 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
//  6 Aug 2001   MGP        Created
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Shell data set: shell identifiers and binding energies

// -------------------------------------------------------------------

#ifndef G4SHELLDATA_HH
#define G4SHELLDATA_HH 1

#include "globals.hh"
#include "g4std/vector"
#include "g4std/map"

class G4DataVector;

class G4ShellData 
{ 
public:

  G4ShellData(G4int minZ = 1, G4int maxZ = 99);

  ~G4ShellData();
 
  size_t NumberOfShells(G4int Z) const;

  G4int ShellId(G4int Z, G4int shellIndex) const;

  const G4DataVector& ShellIdVector(G4int Z) const;

  G4double BindingEnergy(G4int Z, G4int shellIndex) const;

  void LoadData(const G4String& fileName);

  void PrintData() const;

private:

  G4int zMin;
  G4int zMax; 

  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > idMap;
  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > bindingMap;
  G4std::vector<G4int> nShells;

};
 
#endif










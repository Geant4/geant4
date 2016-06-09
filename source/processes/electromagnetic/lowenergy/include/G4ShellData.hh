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
// $Id: G4ShellData.hh,v 1.3 2003/06/16 16:59:48 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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
#include <vector>
#include <map>

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

  // Hide copy constructor and assignment operator 
  G4ShellData& operator=(const G4ShellData& right);
  G4ShellData(const G4ShellData&);

   G4int zMin;
  G4int zMax; 

  std::map<G4int,G4DataVector*,std::less<G4int> > idMap;
  std::map<G4int,G4DataVector*,std::less<G4int> > bindingMap;
  std::vector<G4int> nShells;

};
 
#endif










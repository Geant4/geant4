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

#ifndef G4NuclearLevelStore_hh 
#define G4NuclearLevelStore_hh 1

#include "G4NuclearLevelManager.hh"
#include <map>

class G4NuclearLevelStore
{
private:
  G4NuclearLevelStore();

public:

  static G4NuclearLevelStore* GetInstance();

  G4NuclearLevelManager * GetManager(const G4int Z, const G4int A);


  ~G4NuclearLevelStore();

private:

  G4String GenerateKey(const G4int Z, const G4int A);


  static std::map<G4String,G4NuclearLevelManager*> theManagers;
  static G4String dirName;

};
#endif

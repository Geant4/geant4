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
//
// R&D: Vladimir.Grichine@cern.ch

#ifndef G4VHadDataManager_HH
#define G4VHadDataManager_HH 1

#include "globals.hh"
#include "G4DataVector.hh"

class G4VDataSetAlgorithm;

class G4VHadDataManager 
{ 
public:

  G4VHadDataManager() { };

  virtual ~G4VHadDataManager() { };
 
  virtual G4double FindValue(G4double e, G4int id = 0) const = 0;

  
protected:

private:

};
 
#endif
 











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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef G4NeutronHPDataUsed_h
#define G4NeutronHPDataUsed_h 1

#include "globals.hh"

class G4NeutronHPDataUsed
{
  public:
  
  G4NeutronHPDataUsed()
  {
    theName = "";
    theA = 0;
    theZ = 0;
  }
  
  void SetA(G4double anA){theA = anA;}
  void SetZ(G4int aZ){theZ = aZ;}
  void SetName(G4String aName){theName = aName;}

  G4int GetZ() {return theZ;}
  G4double GetA() {return theA;}
  G4String GetName() {return theName;}
  
  private:
  
  G4String theName;
  G4double theA;
  G4int theZ;
};

#endif

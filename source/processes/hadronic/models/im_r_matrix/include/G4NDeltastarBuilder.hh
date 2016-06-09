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
#ifndef G4NDeltastarBuilder_h
#define G4NDeltastarBuilder_h

#include "G4VXResonanceTable.hh"
#include "G4XNDeltastarTable.hh"

class G4NDeltastarBuilder : public G4VXResonanceTable
{
  public: 
  G4NDeltastarBuilder(const G4String & aName, G4XNDeltastarTable & aT)
      : theT(aT), theS(aName)
  {
    //G4cout << aName <<G4endl;
  }
  virtual G4PhysicsVector* CrossSectionTable() const
  {
    return const_cast<G4PhysicsVector*>(theT.CrossSectionTable(theS));
  }
  
  private:
  
  G4XNDeltastarTable & theT;
  G4String theS;
};
#endif

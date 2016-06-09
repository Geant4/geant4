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

#ifndef G4VNuclearField_h
#define  G4VNuclearField_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4V3DNucleus.hh"
#include "G4HadronicException.hh"

class G4VNuclearField
{
public:
  G4VNuclearField(G4V3DNucleus * aNucleus = 0);
  virtual ~G4VNuclearField();

  void SetNucleus(G4V3DNucleus * aNucleus);
  virtual G4double GetField(const G4ThreeVector & aPosition) = 0;
  virtual G4double GetBarrier() = 0;
  virtual G4double GetCoeff() { return 0; }

protected:
  G4V3DNucleus * theNucleus;
  const G4double radius;

private:

  G4VNuclearField(const  G4VNuclearField &right);
  const G4VNuclearField & operator=(const G4VNuclearField & right);
  G4int operator==(const G4VNuclearField & right) const;
  G4int operator!=(const G4VNuclearField & right) const;

};



inline void G4VNuclearField::SetNucleus(G4V3DNucleus * aNucleus)
{
  theNucleus = aNucleus;
}


#endif


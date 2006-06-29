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


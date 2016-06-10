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
// $Id: G3Pos.hh 67982 2013-03-13 10:36:03Z gcosmo $
//
// ----------------------
// Class description:
//
// G3 volume position

// ----------------------
//
// modified by I.Hrivnacova, 13.10.99

#ifndef G3POS_HH
#define G3POS_HH 1

#include "G4ThreeVector.hh"

class G3Pos
{

public:  // with description

  G3Pos(){;}

  G3Pos(G4String M, G4int C, G4ThreeVector* T, G4int R, G4String O);

  G4bool operator == (const G3Pos& g3p) const;

  virtual ~G3Pos();

  G4String& GetMotherName();

  G4int GetIrot();

  G4ThreeVector* GetPos();

  G4int GetCopy();

  G4String& GetOnly();

private:

  G4String _MotherName;   
  G4int _Copy;
  G4ThreeVector* _Position;
  G4int _Irot;
  G4String _Only;
};

#endif

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
// $Id: G3Pos.hh,v 1.10 2001-07-11 09:58:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

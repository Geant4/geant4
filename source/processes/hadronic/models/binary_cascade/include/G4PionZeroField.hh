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
// -------------------------------------------------------------------
//      GEANT 4 class header file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4PionZeroField.hh
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------

#ifndef G4PionZeroField_h
#define  G4PionZeroField_h 1

#include "G4VNuclearField.hh"
#include "G4V3DNucleus.hh"

class G4PionZeroField: public G4VNuclearField
{
public:
  G4PionZeroField(G4V3DNucleus * nucleus, G4double coeff = 0.042*fermi);
  virtual ~G4PionZeroField();

private:
  G4PionZeroField(const  G4PionZeroField &right);
  const G4PionZeroField & operator=(const G4PionZeroField & right);
  int operator==(const G4PionZeroField & right) const;
  int operator!=(const G4PionZeroField & right) const;

public:
  virtual G4double GetField(const G4ThreeVector & aPosition);
  virtual G4double GetBarrier();
  virtual G4double GetCoeff() { return theCoeff; }

private:
  G4double theCoeff;
};

#endif



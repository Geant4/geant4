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
//      File name:     G4AntiProtonField.hh
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------

#ifndef G4AntiProtonField_h
#define  G4AntiProtonField_h 1

//#include "globals.hh"
#include "G4VNuclearField.hh"
#include "G4V3DNucleus.hh"

class G4AntiProtonField: public G4VNuclearField
{
public:
  G4AntiProtonField(G4V3DNucleus * nucleus, G4double coeff = 1.53*fermi);
  virtual ~G4AntiProtonField();

private:
  G4AntiProtonField(const  G4AntiProtonField &right);
  const G4AntiProtonField & operator=(const G4AntiProtonField & right);
  int operator==(const G4AntiProtonField & right) const;
  int operator!=(const G4AntiProtonField & right) const;

public:
  virtual G4double GetField(const G4ThreeVector & aPosition);
  virtual G4double GetBarrier();
  virtual G4double GetCoeff() { return theCoeff; }

private:
  G4double theCoeff;
};

#endif




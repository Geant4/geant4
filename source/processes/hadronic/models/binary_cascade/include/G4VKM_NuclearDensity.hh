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
//      File name:     G4VKM_NuclearDensity.hh
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------

#ifndef G4VKM_NuclearDensity_hh
#define G4VKM_NuclearDensity_hh

#include "globals.hh"
#include "G4ThreeVector.hh"

class G4VKM_NuclearDensity
{

public:
  G4VKM_NuclearDensity();
  virtual ~G4VKM_NuclearDensity();

  virtual G4double GetDensity(const G4ThreeVector & point) = 0;
  virtual G4double GetDeriv(const G4ThreeVector & point) = 0;
};

inline G4VKM_NuclearDensity::G4VKM_NuclearDensity()
{ }

inline G4VKM_NuclearDensity::~G4VKM_NuclearDensity()
{ }

#endif

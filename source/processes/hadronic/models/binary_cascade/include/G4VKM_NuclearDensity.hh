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

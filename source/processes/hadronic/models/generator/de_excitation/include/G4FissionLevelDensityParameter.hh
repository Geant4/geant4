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
//
// $Id: G4FissionLevelDensityParameter.hh,v 1.6 2001/08/01 17:04:28 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//



#ifndef G4FissionLevelDensityParameter_h
#define G4FissionLevelDensityParameter_h 1


#include "G4VLevelDensityParameter.hh"
#include "G4EvaporationLevelDensityParameter.hh"


class G4FissionLevelDensityParameter : public G4VLevelDensityParameter
{
public:
  G4FissionLevelDensityParameter() {};
  virtual ~G4FissionLevelDensityParameter() {};

private:  
  G4FissionLevelDensityParameter(const G4FissionLevelDensityParameter &right);

  const G4FissionLevelDensityParameter & operator=(const G4FissionLevelDensityParameter &right);
  G4bool operator==(const G4FissionLevelDensityParameter &right) const;
  G4bool operator!=(const G4FissionLevelDensityParameter &right) const;
  
public:
  G4double LevelDensityParameter(const G4int A,const G4int Z,const G4double U) const;

  
private:
  
  G4EvaporationLevelDensityParameter theEvaporationLevelDensityParameter;

};


#endif

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
// $Id: G4VLevelDensityParameter.hh,v 1.6 2002/12/12 19:17:14 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//



#ifndef G4VLevelDensityParameter_h
#define G4VLevelDensityParameter_h 1


#include "globals.hh"

class G4VLevelDensityParameter 
{
public:
  G4VLevelDensityParameter() {};
  virtual ~G4VLevelDensityParameter() {};

private:  
  G4VLevelDensityParameter(const G4VLevelDensityParameter &right);

  const G4VLevelDensityParameter & operator=(const G4VLevelDensityParameter &right);
  G4bool operator==(const G4VLevelDensityParameter &right) const;
  G4bool operator!=(const G4VLevelDensityParameter &right) const;
  
public:
  virtual G4double LevelDensityParameter(const G4int A,const G4int Z,const G4double U) const = 0;

};


#endif

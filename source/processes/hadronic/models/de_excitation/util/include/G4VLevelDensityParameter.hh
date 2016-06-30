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
// $Id: G4VLevelDensityParameter.hh 96634 2016-04-27 09:31:49Z gcosmo $
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

  explicit G4VLevelDensityParameter() {};
  virtual ~G4VLevelDensityParameter() {};

  virtual G4double 
  LevelDensityParameter(G4int A, G4int Z, G4double U) const = 0;

private:  
  G4VLevelDensityParameter(const G4VLevelDensityParameter &right) = delete;
  const G4VLevelDensityParameter & operator=(const G4VLevelDensityParameter &right) = delete;
  G4bool operator==(const G4VLevelDensityParameter &right) const = delete;
  G4bool operator!=(const G4VLevelDensityParameter &right) const = delete;
  
};


#endif

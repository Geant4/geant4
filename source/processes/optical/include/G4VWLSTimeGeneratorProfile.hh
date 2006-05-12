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
//
// GEANT4 Class file
//
//
// File name:     G4VWLSTimeGeneratorProfile.hh
//
// Author:        Pedro Rodrigues, Andreia Trindade
//            
// 
// Creation date: 2006-05-07
//
// Modifications: 
//
// Class Description: 
//
// Abstract class for a WLSTimeGeneratorProfile

// -------------------------------------------------------------------
//

#ifndef G4VWLSTimeGeneratorProfile_h
#define G4VWLSTimeGeneratorProfile_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4MaterialPropertiesTable.hh"

class G4VWLSTimeGeneratorProfile 
{

public:

  G4VWLSTimeGeneratorProfile(const G4String& name);

  virtual ~G4VWLSTimeGeneratorProfile();

  virtual G4double GenerateTime(const G4double time_constant) = 0;
  virtual G4double GenerateTime(const G4MaterialPropertiesTable*) = 0;

protected:

private:

  // hide assignment operator

     G4VWLSTimeGeneratorProfile & operator=
                         (const  G4VWLSTimeGeneratorProfile &right);
     G4VWLSTimeGeneratorProfile(const  G4VWLSTimeGeneratorProfile&);

};

#endif

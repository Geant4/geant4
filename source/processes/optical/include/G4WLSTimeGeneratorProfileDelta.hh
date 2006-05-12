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
// File name:     G4WLSTimeGeneratorProfileDelta.hh
//
// Author:        Pedro Rodrigues, Andreia Trindade
// 
// Creation date: 2006-05-07
//
// Modifications: 
//
// Class Description: Discrete Class of WLSTimeGeneratorProfile
//
// -------------------------------------------------------------------
//

#ifndef G4WLSTimeGeneratorProfileDelta_h
#define G4WLSTimeGeneratorProfileDelta_h 1

#include "G4VWLSTimeGeneratorProfile.hh"

class G4WLSTimeGeneratorProfileDelta : public G4VWLSTimeGeneratorProfile
{

public:

  G4WLSTimeGeneratorProfileDelta(const G4String& name);

  ~G4WLSTimeGeneratorProfileDelta();

  G4double GenerateTime(const G4double time_constant);

  G4double GenerateTime(const G4MaterialPropertiesTable*);

protected:

private:

  // hide assignment operator
 
     G4WLSTimeGeneratorProfileDelta & operator=
                             (const  G4WLSTimeGeneratorProfileDelta &right);
     G4WLSTimeGeneratorProfileDelta(const  G4WLSTimeGeneratorProfileDelta&);

};

#endif

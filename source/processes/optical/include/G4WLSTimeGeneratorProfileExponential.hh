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
// File name:     G4WLSTimeGeneratorProfileExponential.hh
//
// Author:        Pedro Rodrigues, Andreia Trindade
// 
// Creation date: 2006-05-07
//
// Modifications:  
//
// Class Description: 
//

// -------------------------------------------------------------------
//

#ifndef G4WLSTimeGeneratorProfileExponential_h
#define G4WLSTimeGeneratorProfileExponential_h 1

#include "G4VWLSTimeGeneratorProfile.hh"

class G4WLSTimeGeneratorProfileExponential : public G4VWLSTimeGeneratorProfile
{

public:

  G4WLSTimeGeneratorProfileExponential(const G4String& name);

  ~G4WLSTimeGeneratorProfileExponential();

  G4double GenerateTime(const G4double time_constant);

  G4double GenerateTime(const G4MaterialPropertiesTable*);

protected:

private:

  // hide assignment operator
 
     G4WLSTimeGeneratorProfileExponential & operator=
                        (const  G4WLSTimeGeneratorProfileExponential &right);
     G4WLSTimeGeneratorProfileExponential(const  G4WLSTimeGeneratorProfileExponential&);

};

#endif

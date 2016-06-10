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
// $Id: G4tgbPlaceParamLinear.hh 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgbPlaceParamLinear
//
// Class description:
//
// Class to represent simple Cartesian parameterisations along a line.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgbPlaceParamLinear_H
#define G4tgbPlaceParamLinear_H 1

#include "globals.hh"

#include "G4tgbPlaceParameterisation.hh"
#include "G4ThreeVector.hh"

class G4VPhysicalVolume;
class G4tgrPlaceParameterisation;

class G4tgbPlaceParamLinear : public G4tgbPlaceParameterisation
{ 
  public:  // with description

    G4tgbPlaceParamLinear( G4tgrPlaceParameterisation* );
    ~G4tgbPlaceParamLinear();
    void ComputeTransformation(const G4int copyNo,
                               G4VPhysicalVolume *physVol) const;
  private:

    G4ThreeVector theDirection;
    G4double theStep;
    G4double theOffset;
};

#endif

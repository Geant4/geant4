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
// $Id: G4tgbPlaceParamCircle.hh 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgbPlaceParamCircle
//
// Class description:
//
// Class to represent simple Cartesian parameterisations along a circle.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgbPlaceParamCircle_H
#define G4tgbPlaceParamCircle_H 1

#include "globals.hh"

#include "G4VPVParameterisation.hh"
#include "G4tgbPlaceParameterisation.hh"
#include "G4ThreeVector.hh"

class G4tgrPlaceParameterisation;
class G4VPhysicalVolume;

class G4tgbPlaceParamCircle : public G4tgbPlaceParameterisation
{ 
  public:  // with description

    G4tgbPlaceParamCircle( G4tgrPlaceParameterisation* );
    ~G4tgbPlaceParamCircle();

    void ComputeTransformation(const G4int copyNo,
                               G4VPhysicalVolume *physVol) const;

  private:

    void GetNormalToAxis();

  private:

    G4double theRadius;
    G4ThreeVector theCircleAxis;
    G4ThreeVector theDirInPlane;

    G4double theStep;
    G4double theOffset;
};

#endif

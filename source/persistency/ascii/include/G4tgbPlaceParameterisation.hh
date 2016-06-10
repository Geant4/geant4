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
// $Id: G4tgbPlaceParameterisation.hh 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgbPlaceParamSquare
//
// Class description:
//
// Class to represent parameterisations of placements.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgbPlaceParameterisation_H
#define G4tgbPlaceParameterisation_H 1

#include "globals.hh"
#include "geomdefs.hh"
#include "G4VPVParameterisation.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4tgrUtils.hh"

class G4VPhysicalVolume;
class G4tgrPlaceParameterisation;

class G4tgbPlaceParameterisation : public G4VPVParameterisation
{ 
  public:  // with description

    G4tgbPlaceParameterisation( G4tgrPlaceParameterisation* tgrParam );
    virtual ~G4tgbPlaceParameterisation();

    virtual void ComputeTransformation(const G4int, G4VPhysicalVolume *) const;

    void CheckNExtraData( G4tgrPlaceParameterisation* tgrParam,
                          G4int nWcheck, WLSIZEtype st,
                          const G4String& methodName );

    G4int GetNCopies() const { return theNCopies; }
    EAxis GetAxis() const { return theAxis; }

  protected:

    G4int theNCopies;
    EAxis theAxis;
    G4ThreeVector theTranslation;
    G4RotationMatrix* theRotationMatrix;

};

#endif

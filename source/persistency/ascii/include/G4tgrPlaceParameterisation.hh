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
// $Id: G4tgrPlaceParameterisation.hh 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrPlaceParameterisation
//
// Class description:
//
// Class to descripe the positioning of a G4tgrVolume inside another
// G4tgrVolume as a parameterised volume. Several types are possible:
// - Parameterisation of the position and rotation for each copy
// - Parameterisation also of the dimensions
// - Parameterisation of the solid type
// Data is just stored in this class, without any calculation of the
// positions of each copy 
// :POS_PARAM "volu_name" copyNo "parent_name" "parametrisation_type"
//            number_copies step offset extra_data(n words).

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrPlaceParameterisation_h
#define G4tgrPlaceParameterisation_h

#include "globals.hh"

#include <vector>

#include "G4tgrPlace.hh"

class G4tgrPlaceParameterisation : public G4tgrPlace
{
  public:  // with description

    G4tgrPlaceParameterisation();
   ~G4tgrPlaceParameterisation();

    G4tgrPlaceParameterisation( const std::vector<G4String>& p);
      // Creates an object passing the parameters

    // Access functions

    const G4String& GetParamType() const { return theParamType; }
     // GetType returns placement type
    std::vector<G4double> GetExtraData() const { return theExtraData; }
    const G4String& GetRotMatName() const { return theRotMatName; }

    friend std::ostream& operator<<(std::ostream& os,
                                    const G4tgrPlaceParameterisation& obj);
  private:

    G4String theParamType;
    std::vector<G4double> theExtraData;
      // Extra data not common to all parameterisations

    G4String theRotMatName;
      // The rotation matrix (by name, as the rotations
      // matrices are not yet created)
};

#endif

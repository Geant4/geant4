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
// $Id: G4tgrVolumeDivision.hh 68052 2013-03-13 14:38:53Z gcosmo $
//
//
// class G4tgrVolumeDivision
//
// Class description:
//
// Class to manage the geometry info of detector unit made by dividing
// a mother detector unit.
// Construct a volume with the data read from a ':DIV_...' tag:
//   :DIV_NUM  "volu_name" "parent_name" "material" number_divisions
//             "replica_axis_name" offset
//   :DIV_STEP "volu_name" "parent_name" "material" step
//             "replica_axis_name" offset

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrVolumeDivision_h
#define G4tgrVolumeDivision_h

#include "globals.hh"

#include <map>

#include "G4tgrVolume.hh"
#include "G4tgrPlaceDivRep.hh"

typedef std::multimap< G4String, G4String > G4mmss;

//---------------------------------------------------------------------------- 
class G4tgrVolumeDivision : public G4tgrVolume 
{
  public:  // with description

    G4tgrVolumeDivision( const std::vector<G4String>& wl );
   ~G4tgrVolumeDivision();

    // G4bool SetSolid(G4tgrVolume* parentDU, G4bool byStep,
    //                EAxis axis, G4double div_step, G4double offset);
      // Set the solid type and parameters dividing the mother volune

    G4tgrPlaceDivRep* GetPlaceDivision() { return thePlaceDiv; }

    friend std::ostream& operator<<(std::ostream& os,
                                    const G4tgrVolumeDivision& obj);
  private:

    G4tgrPlaceDivRep* thePlaceDiv;
};

#endif

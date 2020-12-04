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
// G4tgrVolumeDivision
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

// Author: P.Arce, CIEMAT (November 2007)
// --------------------------------------------------------------------
#ifndef G4tgrVolumeDivision_hh
#define G4tgrVolumeDivision_hh 1

#include <map>

#include "globals.hh"
#include "G4tgrVolume.hh"
#include "G4tgrPlaceDivRep.hh"

using G4mmss = std::multimap<G4String, G4String>;

class G4tgrVolumeDivision : public G4tgrVolume
{
  public:

    G4tgrVolumeDivision(const std::vector<G4String>& wl);
    ~G4tgrVolumeDivision();

    G4tgrPlaceDivRep* GetPlaceDivision() { return thePlaceDiv; }

    friend std::ostream& operator<<(std::ostream& os,
                                    const G4tgrVolumeDivision& obj);

  private:

    G4tgrPlaceDivRep* thePlaceDiv = nullptr;
};

#endif

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
// $Id: G4tgbVolume.hh 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgbVolume
//
// Class description:
//
// Class to manage the geometry info of any detector unit. The detector units
// created in this class are essentially transient copies of Geant4 physical
// volumes. Thus, they are characterized by a name and the parameters of a
// Geant4 physical volume. 
// They have associated several detector positions, that can be instances of
// G4tgrPlace, G4tgrPlaceDivRep or G4tgrPlaceParameterisation.
// Each detector positioning is done inside a parent. As there can be several
// parents, we will write one parent for each detector position, even if that
// means that parents are repeated.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgbVolume_h
#define G4tgbVolume_h

#include "globals.hh"

#include <vector>
#include <map>

#include "G4tgrVolume.hh"
#include "geomdefs.hh"

class G4tgrPlace;
class G4tgrSolid;
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4AssemblyVolume;

class G4tgbVolume
{
  public:  // with description

    G4tgbVolume();
   ~G4tgbVolume();
    G4tgbVolume( G4tgrVolume* vol);

    void ConstructG4Volumes( const G4tgrPlace* place,
                             const G4LogicalVolume* parentLV );
      // Construct the G4VSolid, G4LogicalVolume and the G4VPhysicalVolume
      // of copy 'copyNo'

    G4VSolid* FindOrConstructG4Solid( const G4tgrSolid* vol);
      // Construct the G4VSolid from the data of the corresponding G4tgrVolume. 
      // Allow to use data from another G4tgrVolume, needed by Boolean solids
      // (that have to construct two solids and then do the Boolean operation)

    G4LogicalVolume* ConstructG4LogVol( const G4VSolid* solid );
      // Construct the G4LogicalVolume and then call the construction of
      // volumes that are positioned inside this LV

    G4VPhysicalVolume* ConstructG4PhysVol( const G4tgrPlace* place,
                                           const G4LogicalVolume* currentLV,
                                           const G4LogicalVolume* parentLV );
      // Construct the G4VPhysicalVolume placing 'curentLV' with position
      // given by the G4tgrPlace 'copyNo' inside 'parentLV'

    void SetCutsInRange( G4LogicalVolume* logvol,
                         std::map<G4String,G4double> cuts );
    void SetCutsInEnergy( G4LogicalVolume* logvol,
                         std::map<G4String,G4double> cuts );


    void CheckNoSolidParams( const G4String& solidType,
                             const unsigned int NoParamExpected,
                             const unsigned int NoParam );
      // Before building a solid of type 'solydType', check if the number
      // of paramenters is the expected one

  G4VSolid* BuildSolidForDivision( G4VSolid* parentSolid, EAxis axis );

    const G4String& GetName() const { return theTgrVolume->GetName(); }
    G4bool GetVisibility() const { return theTgrVolume->GetVisibility(); }
    const G4double* GetColour() const { return theTgrVolume->GetColour(); }

  private:
 
    G4tgrVolume* theTgrVolume;
      // The G4tgrVolume to which it corresponds

    G4AssemblyVolume* theG4AssemblyVolume;
};

#endif 

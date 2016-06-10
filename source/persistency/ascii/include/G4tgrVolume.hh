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
// $Id: G4tgrVolume.hh 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrVolume
//
// Class description:
//
// Abstract base class to manage the geometry info of any volume.
// Volumes created in this class contain the information of a detector volume.
// They have associated several detector placements that can be instances of
// G4tgrPlace, G4tgrPlaceDivision, G4tgrPlaceDivRep or
// G4tgrPlaceParameterisation.
// Each detector positioning is done inside a parent. As there can be several
// parents, one parent for each volume placement will be written, even if that
// means that parents are repeated...

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrVolume_h
#define G4tgrVolume_h

#include "globals.hh"

#include <vector>
#include <map>

class G4tgrSolid;
class G4tgrPlace;
class G4tgrPlaceDivRep;
class G4tgrPlaceParameterisation;

class G4tgrVolume
{
  public:  // with description

    G4tgrVolume();
    G4tgrVolume( const std::vector<G4String>& wl );
    G4tgrVolume( const G4tgrVolume& vol );
    virtual ~G4tgrVolume();

    virtual G4tgrPlace* AddPlace( const std::vector<G4String>& wl );
      // Add a position with the data read from a ':place' tag

    G4tgrPlaceDivRep* AddPlaceReplica( const std::vector<G4String>& wl );
      // Add a replicated position

    G4tgrPlaceParameterisation* AddPlaceParam(const std::vector<G4String>& wl);
      // Add a parameterised position

    void AddVisibility( const std::vector<G4String>& wl );
      // Add visibility flag

    void AddRGBColour( const std::vector<G4String>& wl );
      // Add colour

    void AddCheckOverlaps( const std::vector<G4String>& wl );
      // Add check overlaps flag

    // Accessors

    const G4String& GetName() const {return theName;}
    void SetName(const G4String& name) {theName = name;}
    const G4String& GetType() const {return theType;}
    G4tgrSolid* GetSolid() const {return theSolid;}
    const G4String& GetMaterialName() const {return theMaterialName;}

    const std::vector<G4tgrPlace*> GetPlacements() const {return thePlacements;}
    G4bool GetVisibility() const {return theVisibility;}
    G4double* GetColour() const {return theRGBColour;}
    G4double* GetRGBColour() const {return theRGBColour;}

    G4bool GetCheckOverlaps() const {return theCheckOverlaps;}

    virtual G4tgrVolume* GetVolume( G4int ii ) const;

    friend std::ostream& operator<<(std::ostream& os, const G4tgrVolume& obj);

  protected:   

    G4String theName;   
      // Name of the volume
    G4String theType;   
      // Type of the volume    
    G4String theMaterialName;   
      // Material of which the corresponding PV will be made of
    G4tgrSolid* theSolid;
      // Solid 
    std::vector<G4tgrPlace*> thePlacements;
      // Vector of placements 

    G4bool theVisibility;
    G4double* theRGBColour;
    G4bool theCheckOverlaps;
};

#endif

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
// $Id: G4PhysicalTouchable.hh,v 1.2 2005-02-16 18:14:31 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4PhysicalTouchable
//
// Class description:
//
// This is a concrete class that combines a simple physical volume,  
// with a touchable for the parent volume
//
// Behaviour expected: 
//    - pass all 'VPhysicalVolume' methods to pCurrentVol
//    - respond to GetParentTouchable method
//
// History:
// 14.02.05 J.Apostolakis Created

#ifndef G4PHYSICALTOUCHABLE_HH
#define G4PHYSICALTOUCHABLE_HH

#include "G4Types.hh"
#include "G4String.hh"

#include "geomdefs.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"

class G4LogicalVolume;
class G4VPVParameterisation;
class G4VTouchable;

class G4PhysicalTouchable : public G4VPhysicalVolume
{
  public:  // with description

    G4PhysicalTouchable( G4VPhysicalVolume* pCurrentVol, 
			 const G4VTouchable* pParentTouchable); 

    virtual ~G4PhysicalTouchable();
      // Destructor, will be subclassed. Removes volume from volume Store.

    const G4VTouchable* GetTouchable() const; 
      // Provide parent touchable
    G4VPhysicalVolume* GetCurrentVolume(); 
    const G4VPhysicalVolume* GetCurrentVolume() const; 
      // Access 

    void SetCurrentVolume( G4VPhysicalVolume* pCurrentVol ); 
      // Revise current volume pointer

    inline G4VPhysicalVolume* operator ->(); 
      // 

  private:
    G4VPhysicalVolume*   fpPhysVol; 
      // Current volume pointer

    const G4VTouchable*  fpTouchable;

}; 

// Inline methods

inline G4VPhysicalVolume* G4PhysicalTouchable::operator ->() // const
{
    return( fpPhysVol ? fpPhysVol : 0 );
}
  
#endif

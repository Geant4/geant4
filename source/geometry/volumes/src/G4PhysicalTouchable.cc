// $Id: G4PhysicalTouchable.cc,v 1.1 2005-02-23 10:53:53 japost Exp $
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

#include "G4VPhysicalVolume.hh"
#include "G4PhysicalTouchable.hh"

G4PhysicalTouchable::G4PhysicalTouchable( G4VPhysicalVolume* pCurrentVol, const G4VTouchable* pParentTouchable)

  : G4VPhysicalVolume( 0,                              // G4RotationMatrix *pRot,
		       pCurrentVol->GetTranslation(),  // const G4ThreeVector &tlate,
		       pCurrentVol->GetName(),         // const G4String &pName,
                       pCurrentVol->GetLogicalVolume(), // Current Log Vol
                       0 )                              // physical Mother
{
   if( fpPhysVol ) fpPhysVol = pCurrentVol; 
   fpTouchable= pParentTouchable; 
}

G4PhysicalTouchable::~G4PhysicalTouchable()
      // Destructor, will be subclassed. Removes volume from volume Store.
{
   // ???
}

const G4VTouchable* G4PhysicalTouchable::GetTouchable() const
      // Provide touchable
{
   return fpTouchable; 
}

void G4PhysicalTouchable::SetCurrentVolume( G4VPhysicalVolume* pCurrentVol )
      // Revise current volume pointer
{
   if( fpPhysVol ) fpPhysVol = pCurrentVol; 
   // else{ G4Exception( "Cannot set PhysicalTouchable's current pointer to 0."); 
}

// $Id: G4PhysicalTouchable.cc,v 1.2 2005-03-03 17:08:53 japost Exp $
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
#include "G4TouchableHistory.hh"

G4PhysicalTouchable::G4PhysicalTouchable( G4VPhysicalVolume* pCurrentVol, const G4VTouchable* pParentTouchable)

  : G4VPhysicalVolume( 0,    // G4RotationMatrix *pRot
		       G4ThreeVector(0.,0.,0.), // const G4ThreeVector &tlate,
		       G4String("Not Set"), 
                       0,                       // Current Log Vol
                       0 )                      // physical Mother
{
   fpPhysVol = pCurrentVol; 

   // Set this volume's attributes from the original one.
   if( pCurrentVol ) {
     // Protect copying of attributes - do not use the Physical Volume constructor
     CopyAttributes( pCurrentVol ); 
   }
   else
   {
     G4Exception( "G4PhysicalTouchable constructor( physVol*, const vTouchable*",
		  "InvalidPhysicalVolumePointer", FatalException,
		  "Cannot set PhysicalTouchable's current pointer to Null"); 
   }

   fpTouchable= pParentTouchable; 
   fCreatedParentTouch= false; 
}

G4PhysicalTouchable::G4PhysicalTouchable( G4VPhysicalVolume* pCurrentVol, 
					  const G4NavigationHistory& parentHist)

  : G4VPhysicalVolume( 0,    // G4RotationMatrix *pRot
		       G4ThreeVector(0.,0.,0.), // const G4ThreeVector &tlate,
		       G4String("Not Set"), 
                       0,                       // Current Log Vol
                       0 )                      // physical Mother
{
   fpPhysVol = pCurrentVol; 

   fCreatedParentTouch= true; 
   fpTouchable= new G4TouchableHistory( parentHist ); 

   // Set this volume's attributes from the original one.
   if( pCurrentVol ) {
     // Do it here instead of in the Physical Volume constructor
     CopyAttributes( pCurrentVol ); 
   }
   else
   {
     G4Exception( "G4PhysicalTouchable constructor( physVol*, const navHist*",
		  "InvalidPhysicalVolumePointer", FatalException,
		  "Cannot set PhysicalTouchable's current pointer to Null"); 
   }

}

// Not to be used 
G4PhysicalTouchable::G4PhysicalTouchable
( const G4PhysicalTouchable & )
  : G4VPhysicalVolume( 0,    // G4RotationMatrix *pRot
		       G4ThreeVector(0.,0.,0.), // const G4ThreeVector &tlate,
		       G4String("Not Set"), 
                       0,                       // Current Log Vol
                       0 )                      // physical Mother
{
  G4Exception("G4PhysicalTouchable::G4PhysicalTouchable()",
                  "InvalidCopyConstructor", FatalException,
                  "Copy Constructor is private and not implemented !" ); 
}

void
G4PhysicalTouchable::CopyAttributes( G4VPhysicalVolume* pCurrentVol )
{
  if( pCurrentVol ){
    this->SetRotation(      pCurrentVol->GetRotation() ); 
    this->SetTranslation(   pCurrentVol->GetTranslation() ); 
    this->SetName(          pCurrentVol->GetName() );  
    this->SetLogicalVolume( pCurrentVol->GetLogicalVolume()); 
  }
}

G4PhysicalTouchable::~G4PhysicalTouchable()
      // Parent calls destructor removes volume from volume Store.
{
   if( fCreatedParentTouch ) 
     delete fpTouchable; 
   fpTouchable=0; 
   fCreatedParentTouch= false; 
}

void G4PhysicalTouchable::SetCurrentVolume( G4VPhysicalVolume* pCurrentVol )
      // Revise current volume pointer
{
   if( pCurrentVol ){ fpPhysVol = pCurrentVol; }
   else{ 
     G4Exception( "G4PhysicalTouchable::SetCurrentVolume", 
		  "InvalidPhysicalVolumePointer", FatalException,
		  "Cannot set PhysicalTouchable's current pointer to Null"); 
   }
}

void G4PhysicalTouchable::SetParentTouchable( const G4VTouchable* newParentT ) 
{
   if( fCreatedParentTouch ) {
      delete fpTouchable; 
   }
   fpTouchable= newParentT; 
   fCreatedParentTouch= false; 
}

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
// $Id: G4PhysicalTouchable.cc,v 1.3 2005-07-25 10:02:43 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4PhysicalTouchable implementation
//
// ----------------------------------------------------------------------

#include "G4VPhysicalVolume.hh"
#include "G4PhysicalTouchable.hh"
#include "G4TouchableHistory.hh"

G4PhysicalTouchable::
G4PhysicalTouchable( G4VPhysicalVolume* pCurrentVol,
                     const G4VTouchable* pParentTouchable )
  : G4VPhysicalVolume( 0,                       // G4RotationMatrix *pRot
                       G4ThreeVector(0.,0.,0.), // const G4ThreeVector &tlate,
                       G4String("Not Set"), 
                       0,                       // Current Log Vol
                       0 )                      // physical Mother
{
   fpPhysVol = pCurrentVol; 

   // Set this volume's attributes from the original one.
   //
   if( pCurrentVol )
   {
     // Protect copying of attributes.
     // Do not use the Physical Volume constructor
     //
     CopyAttributes( pCurrentVol ); 
   }
   else
   {
     G4Exception("G4PhysicalTouchable::G4PhysicalTouchable(PV*,const Touch*)",
                 "InvalidPhysicalVolumePointer", FatalException,
                 "Cannot set physical-touchable's current pointer to NULL!"); 
   }
   fpTouchable= pParentTouchable; 
   fCreatedParentTouch= false; 
}

G4PhysicalTouchable::
G4PhysicalTouchable( G4VPhysicalVolume* pCurrentVol, 
                     const G4NavigationHistory& parentHist )

  : G4VPhysicalVolume( 0,                       // G4RotationMatrix *pRot
                       G4ThreeVector(0.,0.,0.), // const G4ThreeVector &tlate,
                       G4String("Not Set"), 
                       0,                       // Current Log Vol
                       0 )                      // physical Mother
{
   fpPhysVol = pCurrentVol; 

   fCreatedParentTouch= true; 
   fpTouchable= new G4TouchableHistory( parentHist ); 

   // Set this volume's attributes from the original one.
   //
   if( pCurrentVol )
   {
     // Do it here instead of in the Physical Volume constructor
     //
     CopyAttributes( pCurrentVol ); 
   }
   else
   {
     G4Exception("G4PhysicalTouchable::G4PhysicalTouchable(PV*,const navHist*)",
                 "InvalidPhysicalVolumePointer", FatalException,
                 "Cannot set physical-touchable's current pointer to NULL!"); 
   }

}

// Not to be used 
//
G4PhysicalTouchable::
G4PhysicalTouchable( const G4PhysicalTouchable & )
  : G4VPhysicalVolume( 0,                       // G4RotationMatrix *pRot
                       G4ThreeVector(0.,0.,0.), // const G4ThreeVector &tlate,
                       G4String("Not Set"), 
                       0,                       // Current Log Vol
                       0 )                      // physical Mother
{
  G4Exception("G4PhysicalTouchable::G4PhysicalTouchable()",
              "InvalidCopyConstructor", FatalException,
              "Copy Constructor is private and not implemented !" ); 
}

void G4PhysicalTouchable::CopyAttributes( G4VPhysicalVolume* pCurrentVol )
{
  if( pCurrentVol )
  {
    this->SetRotation     ( pCurrentVol->GetRotation() ); 
    this->SetTranslation  ( pCurrentVol->GetTranslation() ); 
    this->SetName         ( pCurrentVol->GetName() );  
    this->SetLogicalVolume( pCurrentVol->GetLogicalVolume()); 
  }
}

G4PhysicalTouchable::~G4PhysicalTouchable()
{
  // Parent calls destructor removes volume from volume Store
  //
  if( fCreatedParentTouch ) { delete fpTouchable; }
  fpTouchable = 0; 
  fCreatedParentTouch = false; 
}

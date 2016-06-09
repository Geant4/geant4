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
// $Id: G4FieldManager.cc,v 1.13 2003/11/08 04:08:13 japost Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// -------------------------------------------------------------------

#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4MagneticField.hh"
#include "G4ChordFinder.hh"

G4FieldManager::G4FieldManager(G4Field       *detectorField, 
			       G4ChordFinder *pChordFinder, 
			       G4bool        fieldChangesEnergy
			      )
   : fDetectorField(detectorField), 
     fChordFinder(pChordFinder), 
     fAllocatedChordFinder(false),
     fDefault_Delta_One_Step_Value(0.01*mm), 
     fDefault_Delta_Intersection_Val(0.001*mm),
     fEpsilonMinDefault(5.0e-5), 
     fEpsilonMaxDefault(0.001),
     fEpsilonMin( fEpsilonMinDefault ),
     fEpsilonMax( fEpsilonMaxDefault)
{ 
   fDelta_One_Step_Value= fDefault_Delta_One_Step_Value;
   fDelta_Intersection_Val= fDefault_Delta_Intersection_Val;
   if ( detectorField )
     fFieldChangesEnergy= detectorField->DoesFieldChangeEnergy();
   else
     fFieldChangesEnergy= fieldChangesEnergy;
}

G4FieldManager::G4FieldManager(G4MagneticField *detectorField)
   : fDetectorField(detectorField), fAllocatedChordFinder(true),
     fFieldChangesEnergy(false), 
     fDefault_Delta_One_Step_Value(0.01*mm),
     fDefault_Delta_Intersection_Val(0.001*mm),
     fEpsilonMinDefault(5.0e-5), 
     fEpsilonMaxDefault(0.001),
     fEpsilonMin( fEpsilonMinDefault ),
     fEpsilonMax( fEpsilonMaxDefault)
{
   fChordFinder= new G4ChordFinder( detectorField );
   fDelta_One_Step_Value= fDefault_Delta_One_Step_Value;
   fDelta_Intersection_Val= fDefault_Delta_Intersection_Val;
}

void G4FieldManager::ConfigureForTrack( const G4Track * ) 
{
   // Default is to do nothing!
   ;
}

G4FieldManager::~G4FieldManager()
{
   if( fAllocatedChordFinder ){
      delete fChordFinder;
   }
}

void
G4FieldManager::CreateChordFinder(G4MagneticField *detectorMagField)
{
   if ( fAllocatedChordFinder )
      delete fChordFinder;
   fChordFinder= new G4ChordFinder( detectorMagField );
   fAllocatedChordFinder= true;
}

G4bool G4FieldManager::SetDetectorField(G4Field *pDetectorField)
{
   fDetectorField= pDetectorField;

   if ( pDetectorField )
     fFieldChangesEnergy= pDetectorField->DoesFieldChangeEnergy();
   else
     fFieldChangesEnergy= false;   //  No field 

   return false;
}


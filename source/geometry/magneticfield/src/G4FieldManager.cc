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
// $Id: G4FieldManager.cc,v 1.15 2007-12-07 15:34:10 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------

#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4MagneticField.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManagerStore.hh"

G4FieldManager::G4FieldManager(G4Field       *detectorField, 
			       G4ChordFinder *pChordFinder, 
			       G4bool        fieldChangesEnergy
			      )
   : fDetectorField(detectorField), 
     fChordFinder(pChordFinder), 
     fAllocatedChordFinder(false),
     fEpsilonMinDefault(5.0e-5), 
     fEpsilonMaxDefault(0.001),
     fDefault_Delta_One_Step_Value(0.01*mm), 
     fDefault_Delta_Intersection_Val(0.001*mm),
     fEpsilonMin( fEpsilonMinDefault ),
     fEpsilonMax( fEpsilonMaxDefault)
{ 
   fDelta_One_Step_Value= fDefault_Delta_One_Step_Value;
   fDelta_Intersection_Val= fDefault_Delta_Intersection_Val;
   if ( detectorField )
     fFieldChangesEnergy= detectorField->DoesFieldChangeEnergy();
   else
     fFieldChangesEnergy= fieldChangesEnergy;

   // Add to store
   G4FieldManagerStore::Register(this);
}

G4FieldManager::G4FieldManager(G4MagneticField *detectorField)
   : fDetectorField(detectorField), fAllocatedChordFinder(true),
     fEpsilonMinDefault(5.0e-5), 
     fEpsilonMaxDefault(0.001),
     fFieldChangesEnergy(false), 
     fDefault_Delta_One_Step_Value(0.01*mm),
     fDefault_Delta_Intersection_Val(0.001*mm),
     fEpsilonMin( fEpsilonMinDefault ),
     fEpsilonMax( fEpsilonMaxDefault)
{
   fChordFinder= new G4ChordFinder( detectorField );
   fDelta_One_Step_Value= fDefault_Delta_One_Step_Value;
   fDelta_Intersection_Val= fDefault_Delta_Intersection_Val;
   // Add to store
   G4FieldManagerStore::Register(this);
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
   G4FieldManagerStore::DeRegister(this);
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


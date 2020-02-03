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
// G4FieldManager implementation
//
// Author: John Apostolakis, 10.03.97 - design and implementation
// -------------------------------------------------------------------

#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4MagneticField.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManagerStore.hh"
#include "G4SystemOfUnits.hh"

G4double G4FieldManager::fDefault_Delta_One_Step_Value= 0.01 *    millimeter;
G4double G4FieldManager::fDefault_Delta_Intersection_Val= 0.001 * millimeter;

G4FieldManager::G4FieldManager(G4Field* detectorField, 
                               G4ChordFinder* pChordFinder, 
                               G4bool fieldChangesEnergy )
   : fDetectorField(detectorField), 
     fChordFinder(pChordFinder), 
     fDelta_One_Step_Value( fDefault_Delta_One_Step_Value ), 
     fDelta_Intersection_Val( fDefault_Delta_Intersection_Val ),     
     fEpsilonMin( fEpsilonMinDefault ),
     fEpsilonMax( fEpsilonMaxDefault)
{ 
   if ( detectorField != nullptr )
   {
     fFieldChangesEnergy = detectorField->DoesFieldChangeEnergy();
   }
   else
   {
     fFieldChangesEnergy = fieldChangesEnergy;
   }

   // Add to store
   //
   G4FieldManagerStore::Register(this);
}

G4FieldManager::G4FieldManager(G4MagneticField* detectorField)
   : fDetectorField(detectorField), fAllocatedChordFinder(true),
     fDelta_One_Step_Value(   fDefault_Delta_One_Step_Value ), 
     fDelta_Intersection_Val( fDefault_Delta_Intersection_Val ),
     fEpsilonMin( fEpsilonMinDefault ),
     fEpsilonMax( fEpsilonMaxDefault )
{
   fChordFinder = new G4ChordFinder( detectorField );

   // Add to store
   //
   G4FieldManagerStore::Register(this);
}

G4FieldManager* G4FieldManager::Clone() const
{
    G4Field* aField = nullptr;
    G4FieldManager* aFM = nullptr;
    G4ChordFinder* aCF = nullptr;
    try {
        if ( fDetectorField != nullptr )
        {
           aField = fDetectorField->Clone();
        }

        // Create a new field manager, note that we do not set
        // any chordfinder now
        //
        aFM = new G4FieldManager( aField , nullptr , fFieldChangesEnergy );

        // Check if originally we have the fAllocatedChordFinder variable
        // set, in case, call chord constructor
        //
        if ( fAllocatedChordFinder )
        {
            aFM->CreateChordFinder( dynamic_cast<G4MagneticField*>(aField) );
        }
        else
        {
            // Chord was specified by user, should we clone?
            // TODO: For the moment copy pointer, to be understood
            //       if cloning of ChordFinder is needed
            //
            aCF = fChordFinder; /*->Clone*/
            aFM->fChordFinder = aCF;
        }

        // Copy values of other variables

        aFM->fEpsilonMax = fEpsilonMax;
        aFM->fEpsilonMin = fEpsilonMin;
        aFM->fDelta_Intersection_Val = fDelta_Intersection_Val;
        aFM->fDelta_One_Step_Value = fDelta_One_Step_Value;
          // TODO: Should we really add to the store the cloned FM?
          //       Who will use this?
    }
    catch ( ... )
    {
        // Failed creating clone: probably user did not implement Clone method
        // in derived classes?
        // Perform clean-up after ourselves...
        delete aField;
        delete aFM;
        delete aCF;
        throw;
    }
    return aFM;
}

void G4FieldManager::ConfigureForTrack( const G4Track * ) 
{
   // Default is to do nothing!
}

G4FieldManager::~G4FieldManager()
{
   if( fAllocatedChordFinder )
   {
      delete fChordFinder;
   }
   G4FieldManagerStore::DeRegister(this);
}

void
G4FieldManager::CreateChordFinder(G4MagneticField* detectorMagField)
{
   if ( fAllocatedChordFinder )
   { 
      delete fChordFinder;
   }
   fAllocatedChordFinder = false;

   if( detectorMagField != nullptr )
   {
      fChordFinder = new G4ChordFinder( detectorMagField );
      fAllocatedChordFinder = true;
   }
   else
   {
      fChordFinder = nullptr;
   }
}

void G4FieldManager::InitialiseFieldChangesEnergy()
{
   if ( fDetectorField != nullptr )
   {
     fFieldChangesEnergy = fDetectorField->DoesFieldChangeEnergy();
   }
   else
   {
     fFieldChangesEnergy = false;   // No field, no change!
   }
}

G4bool G4FieldManager::SetDetectorField(G4Field* pDetectorField,
                                        G4int failMode )
{
   G4VIntegrationDriver* driver = nullptr;
   G4EquationOfMotion* equation = nullptr;
   // G4bool  compatibleField = false;
   G4bool ableToSet = false;

   fDetectorField = pDetectorField;
   InitialiseFieldChangesEnergy();
   
   // Must 'propagate' the field to the dependent classes
   //
   if( fChordFinder != nullptr )
   {
     failMode= std::max( failMode, 1) ;
       // If a chord finder exists, warn in case of error!
      
     driver = fChordFinder->GetIntegrationDriver();
     if( driver != nullptr )
     {
       equation = driver->GetEquationOfMotion();

       // Should check the compatibility between the
       // field and the equation HERE

       if( equation != nullptr )
       {
          equation->SetFieldObj(pDetectorField);
          ableToSet = true;
       }
     }
   }
   
   if( !ableToSet && (failMode > 0) )
   {
     // If this fails, report the issue !

     G4ExceptionDescription msg;
     msg << "Unable to set the field in the dependent objects of G4FieldManager"
         << G4endl;
     msg << "All the dependent classes must be fully initialised,"
         << "before it is possible to call this method." << G4endl;
     msg << "The problem encountered was the following: " << G4endl;
     if( fChordFinder == nullptr ) { msg << "  No ChordFinder. " ; }
     else if( driver == nullptr )  { msg << "  No Integration Driver set. ";}
     else if( equation == nullptr ) { msg << "  No Equation found. " ; }
     // else if( !compatibleField ) { msg << "  Field not compatible. ";}
     else { msg << "  Can NOT find reason for failure. ";}
     msg << G4endl;
     G4ExceptionSeverity severity = (failMode != 1)
                                  ? FatalException : JustWarning ;
     G4Exception("G4FieldManager::SetDetectorField", "Geometry001",
                 severity, msg);
   }
   return ableToSet;
}

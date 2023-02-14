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
G4bool   G4FieldManager::fVerboseConstruction= false;

G4double G4FieldManager::fMaxAcceptedEpsilon= 0.01; //  Legacy value.  Future value = 0.001
//   Requesting a large epsilon (max) value provides poor accuracy for
//   every integration segment.
//   Problems occur because some methods (including DormandPrince(7)45 the estimation of local
//   error appears to be a substantial underestimate at large epsilon values ( > 0.001 ).
//   So the value for fMaxAcceptedEpsilon is recommended to be 0.001 or below.

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

   if( fVerboseConstruction)
      G4cout << "G4FieldManager/ctor#1 fEpsilon Min/Max:  eps_min = " << fEpsilonMin << " eps_max=" << fEpsilonMax << G4endl;

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

   if( fVerboseConstruction )
      G4cout << "G4FieldManager/ctor#2 fEpsilon Min/Max:  eps_min = " << fEpsilonMin << " eps_max=" << fEpsilonMax << G4endl;
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

    G4cout << "G4FieldManager/clone fEpsilon Min/Max:  eps_min = " << fEpsilonMin << " eps_max=" << fEpsilonMax << G4endl;
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

G4bool G4FieldManager::SetMaximumEpsilonStep( G4double newEpsMax )
{
  G4bool succeeded= false;
  if(    (newEpsMax > 0.0) && ( newEpsMax <= fMaxAcceptedEpsilon)
     &&  (fMinAcceptedEpsilon <= newEpsMax ) ) // (std::fabs(1.0+newEpsMax)>1.0) )
  {
    if(newEpsMax >= fEpsilonMin){
      fEpsilonMax = newEpsMax;
      succeeded = true;
      if (fVerboseConstruction)
      {
        G4cout << "G4FieldManager/SetEpsMax :  eps_max = " << std::setw(10) << fEpsilonMax
               << " ( Note: unchanged eps_min=" << std::setw(10) << fEpsilonMin << " )" << G4endl;
      }
    } else {
      G4ExceptionDescription erm;
      erm << " Call to set eps_max = " << newEpsMax << " . The problem is that"
          << " its value must be at larger or equal to eps_min= " << fEpsilonMin << G4endl;
      erm << " Modifying both to the same value " << newEpsMax << " to ensure consistency."
          << G4endl
          << " To avoid this warning, please set eps_min first, and ensure that "
          << " 0 < eps_min <= eps_max <= " << fMaxAcceptedEpsilon << G4endl;
      
      fEpsilonMax = newEpsMax;
      fEpsilonMin = newEpsMax;
      G4String methodName = G4String("G4FieldManager::")+ G4String(__func__);
      G4Exception(methodName.c_str(), "Geometry003", JustWarning, erm);
    }
  }
  else
  {
    G4ExceptionDescription erm;
    G4String paramName("eps_max");
    ReportBadEpsilonValue(erm, newEpsMax, paramName );
    G4String methodName = G4String("G4FieldManager::")+ G4String(__func__);
    G4Exception(methodName.c_str(), "Geometry001", FatalException, erm);
  }
  return succeeded;
}

// -----------------------------------------------------------------------------

G4bool G4FieldManager::SetMinimumEpsilonStep( G4double newEpsMin )
{
  G4bool succeeded= false;

  if( fMinAcceptedEpsilon <= newEpsMin  &&  newEpsMin <= fMaxAcceptedEpsilon )
  {
    fEpsilonMin = newEpsMin;
    //*********
    succeeded= true;

    if (fVerboseConstruction)
    {
      G4cout << "G4FieldManager/SetEpsMin :  eps_min = "
             << std::setw(10) << fEpsilonMin << G4endl;
    }
    if( fEpsilonMax < fEpsilonMin ){
      // Ensure consistency
      G4ExceptionDescription erm;
      erm << "Setting eps_min = " << newEpsMin
          << " For consistency set eps_max= " << fEpsilonMin
          << " ( Old value = " << fEpsilonMax << " )" << G4endl;
      fEpsilonMax = fEpsilonMin;
      G4String methodName = G4String("G4FieldManager::")+ G4String(__func__);
      G4Exception(methodName.c_str(), "Geometry003", JustWarning, erm);
    }
  }
  else
  {
    G4ExceptionDescription erm;
    G4String paramName("eps_min");
    ReportBadEpsilonValue(erm, newEpsMin, paramName );
    G4String methodName = G4String("G4FieldManager::")+ G4String(__func__);
    G4Exception(methodName.c_str(), "Geometry001", FatalException, erm);
  }
  return succeeded;
}

// -----------------------------------------------------------------------------

G4double G4FieldManager::GetMaxAcceptedEpsilon()
{
   return fMaxAcceptedEpsilon;
}

// -----------------------------------------------------------------------------

G4bool   G4FieldManager::SetMaxAcceptedEpsilon(G4double maxAcceptValue, G4bool softFailure)
// Set value -- within limits
{
  G4bool success= false;
  // Limit for warning and absolute limit chosen from experience in and 
  // investigation of integration with G4DormandPrince745 in HEP-type setups.
  if( maxAcceptValue <= fMaxWarningEpsilon )
  {
    fMaxAcceptedEpsilon= maxAcceptValue;
    success= true;
  }
  else
  {
    G4ExceptionDescription erm;
    G4ExceptionSeverity  severity;

    G4cout << "G4FieldManager::" << __func__ 
           << " Parameters:   fMaxAcceptedEpsilon = " << fMaxAcceptedEpsilon
           << " fMaxFinalEpsilon = " << fMaxFinalEpsilon << G4endl;
    
    if( maxAcceptValue <= fMaxFinalEpsilon )
    {
      success= true;
      fMaxAcceptedEpsilon = maxAcceptValue;
      // Integration is poor, and robustness will likely suffer
      erm << "Proposed value for maximum-accepted-epsilon = " << maxAcceptValue
          << " is larger than the recommended = " << fMaxWarningEpsilon
          << G4endl
          << "This may impact the robustness of integration of tracks in field."
          << G4endl
          << "The request was accepted and the value = "  << fMaxAcceptedEpsilon
          << " , but future releases are expected " << G4endl
          << " to tighten the limit of acceptable values to " 
          << fMaxWarningEpsilon << G4endl << G4endl          
          << "Suggestion: If you need better performance investigate using "
          << "alternative, low-order RK integration methods or " << G4endl
          << " helix-based methods (for pure B-fields) for low(er) energy tracks, "
          << " especially electrons if you need better performance." << G4endl;
      severity= JustWarning;
    }
    else
    {
      fMaxAcceptedEpsilon= fMaxFinalEpsilon;
      erm << " Proposed value for maximum accepted epsilon " << maxAcceptValue
          << " is larger than the top of the range = " << fMaxFinalEpsilon
          << G4endl;
      if( softFailure )
         erm << " Using the latter value instead." << G4endl;
      erm << G4endl;
      erm << " Please adjust to request maxAccepted <= " << fMaxFinalEpsilon
          << G4endl << G4endl;
      if( softFailure == false )
         erm << " NOTE: you can accept the ceiling value and turn this into a " 
             << " warning by using a 2nd argument  " << G4endl
             << " in your call to SetMaxAcceptedEpsilon:  softFailure = true ";
      severity = softFailure ? JustWarning : FatalException; 
      // if( softFailure ) severity= JustWarning;
      // else severity= FatalException;
      success = false;
    }
    G4String methodName = G4String("G4FieldManager::")+ G4String(__func__);
    G4Exception(methodName.c_str(), "Geometry003", severity, erm);    
  }
  return success;
}

// -----------------------------------------------------------------------------

void G4FieldManager::
ReportBadEpsilonValue(G4ExceptionDescription& erm, G4double value, G4String& name) const
{
  erm << "Incorrect proposed value of " << name << " = " << value << G4endl
      << " Its value is outside the permitted range from "
      << fMinAcceptedEpsilon << "  to " <<  fMaxAcceptedEpsilon << G4endl
      << " Clarification: " << G4endl;
  G4long oldPrec = erm.precision();
  if(value < fMinAcceptedEpsilon )
  {
    erm << "  a) The value must be positive and enough larger than the accuracy limit"
        << " of the (G4)double type - ("
        <<  (value < fMinAcceptedEpsilon ? "FAILED" : "OK" ) << ")" << G4endl
        << "     i.e. std::numeric_limits<G4double>::epsilon()= "
        <<  std::numeric_limits<G4double>::epsilon()
        << " to ensure that integration " << G4endl
        << "     could potentially achieve this acccuracy." << G4endl
        << "     Minimum accepted eps_min/max value = " << fMinAcceptedEpsilon << G4endl;
  }
  else if( value > fMaxAcceptedEpsilon)
  {
    erm << "  b) It must be smaller than (or equal) " << std::setw(8)
        << std::setprecision(4) << fMaxAcceptedEpsilon
        << " to ensure robustness of integration - ("
        << (( value < fMaxAcceptedEpsilon) ? "OK" : "FAILED" ) << ")" << G4endl;
  }
  else
  {
    G4bool badRoundoff = (std::fabs(1.0+value) == 1.0);
    erm << "  Unknown ERROR case -- extra check: " << G4endl;
    erm << "  c) as a floating point number (of type G4double) the sum (1+" << name
        << " ) must be > 1 , ("
        <<  (badRoundoff ? "FAILED" : "OK" ) << ")" << G4endl
        << "     Now    1+eps_min          = " << std::setw(20)
        << std::setprecision(17) << (1+value) << G4endl
        << "     and   (1.0+" << name << ") - 1.0 = " << std::setw(20)
        << std::setprecision(9) << (1.0+value)-1.0;
  }
  erm.precision(oldPrec);
}

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
// $Id: G4FieldManager.hh,v 1.14 2003/12/09 14:44:59 japost Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//  
// class G4FieldManager
//
// Class description:
//
// A class to manage (Store) a pointer to the Field subclass that
// describes the field of a detector (magnetic, electric or other).
// Also stores a reference to the chord finder.
//
// The G4FieldManager class exists to allow the user program to specify 
// the electric, magnetic and/or other field(s) of the detector.
// (OR, in the future, of a part of it - planned to be a logical volume).
// It also stores a pointer to the ChordFinder object that can do the
// propagation in this field. All geometrical track "advancement" 
// in the field is handled by this ChordFinder object.
//
// G4FieldManager allows the other classes/object (of the MagneticField 
// & other class categories) to find out whether a detector field object 
// exists and what that object is.
//
// The Chord Finder must be created either by calling CreateChordFinder
// for a Magnetic Field or by the user creating a  a Chord Finder object
// "manually" and setting this pointer.
//
// A default FieldManager is created by the singleton class
// G4NavigatorForTracking and exists before main is called.
// However a new one can be created and given to G4NavigatorForTracking.
//
// Our current design envisions that one Field manager is 
// valid for each region detector.

// History:
// - 10.03.97 John Apostolakis, design and implementation.
// -------------------------------------------------------------------

#ifndef G4FIELDMANAGER_HH
#define G4FIELDMANAGER_HH 1

#include "globals.hh"

class G4Field;
class G4MagneticField;
class G4ChordFinder;
class G4Track;  // Forward reference for parameter configuration

class G4FieldManager
{
  public:  // with description
     G4FieldManager(G4Field       *detectorField=0, 
		    G4ChordFinder *pChordFinder=0, 
		    G4bool       b=true ); // fieldChangesEnergy is taken from field
          // General constructor for any field.
          // -> Must be set with field and chordfinder for use.
     G4FieldManager(G4MagneticField *detectorMagneticField);
          // Creates ChordFinder
          //   - assumes pure magnetic field (so Energy constant)
     virtual ~G4FieldManager();

     G4bool          SetDetectorField(G4Field *detectorField);
     inline const G4Field*  GetDetectorField() const;
     inline G4bool          DoesFieldExist() const;
        // Set, get and check the field object

     void            CreateChordFinder(G4MagneticField *detectorMagField);
     inline void     SetChordFinder(G4ChordFinder *aChordFinder);
     inline G4ChordFinder*  GetChordFinder();
     inline const G4ChordFinder*  GetChordFinder() const;
        // Create, set or get the associated Chord Finder

     virtual void   ConfigureForTrack( const G4Track * ); 
        // Setup the choice of the configurable parameters 
        //    relying on the current track's energy, particle identity, ..
        //  Note: In addition to the values of member variables, 
        //         a user can use this to change the ChordFinder, the field, ...

  public:  // with description
   
     inline G4double GetDeltaIntersection() const;  // virtual ?
       // Accuracy for boundary intersection.

     inline G4double GetDeltaOneStep() const;      // virtual ?
       // Accuracy for one tracking/physics step.

     inline void     SetAccuraciesWithDeltaOneStep(G4double valDeltaOneStep); 
       // Sets both accuracies, maintaining a fixed ratio for accuracties 
       // of volume Intersection and Integration (in One Step) 

     inline void     SetDeltaOneStep(G4double valueD1step); 
      // Set accuracy for integration of one step.   (only)
     inline void     SetDeltaIntersection(G4double valueDintersection); 
      // Set accuracy of  intersection of a volume.  (only)

     inline G4double  GetMinimumEpsilonStep() const;
     inline void      SetMinimumEpsilonStep( G4double newEpsMin );
     // Minimum for Relative accuracy of a Step 

     inline G4double  GetMaximumEpsilonStep() const;
     inline void      SetMaximumEpsilonStep( G4double newEpsMax );
     // Maximum for Relative accuracy of a Step 
 
     inline G4bool   DoesFieldChangeEnergy() const;
     inline void     SetFieldChangesEnergy(G4bool value);
       //  For electric field this should be true
       //  For magnetic field this should be false

  private:

     G4FieldManager(const G4FieldManager&);
     G4FieldManager& operator=(const G4FieldManager&);
       // Private copy constructor and assignment operator.

  private:

     G4Field*        fDetectorField;
     G4ChordFinder*  fChordFinder;

     G4bool          fAllocatedChordFinder; // Did we used "new" to
					    // create fChordFinder ?
     G4bool          fFieldChangesEnergy;

     //  Values for the required accuracies
     //
     G4double  fDelta_One_Step_Value;      //  for one tracking/physics step
     G4double  fDelta_Intersection_Val;    //  for boundary intersection

     G4double  fDefault_Delta_One_Step_Value;   // = 0.25 * mm;
     G4double  fDefault_Delta_Intersection_Val; // = 0.1 * mm;

     //  Values for the small possible relative accuracy of a step
     //  (corresponding to the greatest possible integration accuracy)

     G4double  fEpsilonMinDefault;   // Can be 1.0e-5 to 1.0e-10 ...
     G4double  fEpsilonMaxDefault;   // Can be 1.0e-3 to 1.0e-8 ...
     G4double  fEpsilonMin; 
     G4double  fEpsilonMax;
};

// Our current design and implementation expect that a particular
// geometrical region has a Field manager.
//  By default a Field Manager is created for the world volume, and
//  will be utilised for all volumes unless it is overridden by a 'local'
//  field manager.

// Note also that a region with both electric E and magnetic B field will 
//  have these treated as one field.
// Similarly it could be extended to treat other fields as additional components
//  of a single field type.


// Implementation of inline functions

#include "G4FieldManager.icc"

#endif   /*  G4FIELDMANAGER_HH */

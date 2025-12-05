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
// G4FieldManager
//
// Class description:
//
// A class to manage (Store) a pointer to the Field subclass that
// describes the field of a detector (magnetic, electric or other).
// Also stores a reference to the chord finder.
//
// The G4FieldManager class exists to allow the user program to specify 
// the electric, magnetic and/or other field(s) of the detector.
// 
// A field manager can be set to a logical volume (or to more than one), 
// in order to vary its field from that of the world.  In this manner
// a zero or constant field can override a global field,  a more or 
// less exact version can override the external approximation, lower
// or higher precision for tracking can be specified, a different 
// stepper can be chosen for different volumes, ...
//
// It also stores a pointer to the ChordFinder object that can do the
// propagation in this field. All geometrical track "advancement" 
// in the field is handled by this ChordFinder object.
//
// G4FieldManager allows the other classes/object (of the MagneticField 
// & other class categories) to find out whether a detector field object 
// exists and what that object is.
//
// The Chord Finder must be created either by calling CreateChordFinder
// for a Magnetic Field or by the user creating a Chord Finder object
// "manually" and setting the pointer.
//
// A default FieldManager is created by the singleton class
// G4NavigatorForTracking and exists before main is called.
// However a new one can be created and given to G4NavigatorForTracking.
//
// Our current design envisions that one Field manager is 
// valid for each region detector.
//
// It is expected that a particular geometrical region has a Field manager.
// By default a Field Manager is created for the world volume, and
// will be utilised for all volumes unless it is overridden by a 'local'
// field manager.
// Note also that a region with both electric E and magnetic B field will 
// have these treated as one field.
// Similarly it could be extended to treat other fields as additional
// components of a single field type.

// Author: John Apostolakis (CERN), 10.03.1997 - Design and implementation
// -------------------------------------------------------------------
#ifndef G4FIELDMANAGER_HH
#define G4FIELDMANAGER_HH

#include "G4FieldParameters.hh"
#include "globals.hh"

class G4Field;
class G4MagneticField;
class G4ChordFinder;
class G4Track;  // Forward reference for parameter configuration

/**
 * @brief G4FieldManager is a manager (store) for a pointer to the Field
 * subclass that describes the field of a detector (magnetic, electric or
 * other). It also stores a reference to the chord finder.
 * A field manager can be set to a logical volume (or to more than one), 
 * in order to vary its field from that of the world volume. In this manner
 * a zero or constant field can override a global field, a more or less exact
 * version can override the external approximation, lower or higher precision
 * for tracking can be specified, a different stepper can be chosen for
 * different volumes, etc...
 * The Chord Finder must be created either by calling CreateChordFinder()
 * for a Magnetic Field or by the user creating a Chord Finder object
 * "manually" and setting the pointer.
 * The current design envisions that one Field manager is valid for each
 * detector region. It is expected that a particular geometrical region has
 * a Field manager. By default a Field Manager is created for the world volume,
 * and will be utilised for all volumes unless it is overridden by a 'local'
 * field manager.
 * Note also that a region with both electric E and magnetic B field will 
 * have these treated as one field. Similarly it could be extended to treat
 * other fields as additional components of a single field type.
 */

class G4FieldManager
{
  public:

    /**
     * General Constructor for any field. Must be set with field and chord finder
     * for use.
     *  @param[in] detectorField Pointer to the field.
     *  @param[in] pChordFinder Pointer to the chord finder object.
     *  @param[in] b Flag to indicate if the field changes the energy; it is
     *             taken from the provided field, if specified.
     */
    G4FieldManager(G4Field* detectorField = nullptr, 
                   G4ChordFinder* pChordFinder = nullptr, 
                   G4bool b = true ); // fieldChangesEnergy is taken from field

    /**
     * Constructor creating the chord finder. It assumes pure magnetic field,
     * so energy constant.
     *  @param[in] detectorMagneticField Pointer to the magnetic field.
     */
    G4FieldManager(G4MagneticField* detectorMagneticField);

    /**
     * Virtual Destructor.
     */
    virtual ~G4FieldManager();

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4FieldManager(const G4FieldManager&) = delete;
    G4FieldManager& operator=(const G4FieldManager&) = delete;

    /**
     * Pushes the field to the equation. Failure to push the field (due to
     * absence of a chord finder, driver, stepper or equation) is
     *  - '0' = quiet      : Do not complain if chordFinder == 0
     *                       (It will still warn for other error);
     *  - '1' = warn       : a warning if anything is missing;
     *  - '2'/else = FATAL : a fatal error for all other values.
     *  @param[in] detectorField Pointer to the field.
     *  @param[in] failMode Flag (0/1/2) for selected failure mode.
     *  @returns Success (true) or failure (false).
     */
    G4bool SetDetectorField(G4Field* detectorField, G4int failMode = 0);

    /**
     * Pushes the field to this class only -- no further.
     * Should be used to initialise this field, only *before* creating
     * the chord finder and its dependent classes.
     * User is then responsible to ensure that:
     *    i) an equation, stepper, driver and chord finder are created;
     *   ii) this field is used by the equation.
     *  @param[in] detectorField Pointer to the field.
     */
    inline void ProposeDetectorField(G4Field* detectorField);

    /**
     * Pushes the field to the equation and keeps its address.
     * Can be used only once the equation, stepper, driver and chord finder
     * have all been created; else it is an error.
     *  @param[in] detectorField Pointer to the field.
     */
    inline void ChangeDetectorField(G4Field* detectorField);    
        
    /**
     * Methods to get and check (existance of) the field object.
     */
    inline const G4Field* GetDetectorField() const;
    inline G4bool DoesFieldExist() const;

    /**
     * Methods to create, set or get the associated Chord Finder.
     */
    void CreateChordFinder(G4MagneticField* detectorMagField);
    inline void SetChordFinder(G4ChordFinder* aChordFinder);
    inline G4ChordFinder* GetChordFinder();
    inline const G4ChordFinder* GetChordFinder() const;

    /**
     * Setups the choice of the configurable parameters, relying on the
     * current track's energy, particle identity...
     * Note: in addition to the values of member variables, a user can use
     * this to change the ChordFinder, the field, etc.
     *  @param[in] pTrack Pointer to a track.
     */
    virtual void ConfigureForTrack( const G4Track* pTrack );

    /**
     * Static methods to set/get the global field.
     */
    static void SetGlobalFieldManager(G4FieldManager* fieldManager);
    static G4FieldManager* GetGlobalFieldManager();

    /**
     * Returns the accuracy for boundary intersection.
     */
    inline G4double GetDeltaIntersection() const;

    /**
     * Returns the accuracy for one tracking/physics step.
     */
    inline G4double GetDeltaOneStep() const;

    /**
     * Sets both accuracies, maintaining a fixed ratio for accuracies 
     * of volume Intersection and Integration (in One Step).
     */
    inline void SetAccuraciesWithDeltaOneStep(G4double valDeltaOneStep); 

    /**
     * Sets the accuracy for integration of one step (only).
     */
    inline void SetDeltaOneStep(G4double valueD1step); 

    /**
     * Sets the accuracy of intersection of a volume (only).
     */
    inline void     SetDeltaIntersection(G4double valueDintersection); 

    /**
     * Methods to set/get the minimum for Relative accuracy of a Step.
     */
    inline G4double GetMinimumEpsilonStep() const;
    G4bool SetMinimumEpsilonStep( G4double newEpsMin );

    /**
     * Methods to set/get the maximum for Relative accuracy of a Step.
     */
    inline G4double GetMaximumEpsilonStep() const;
    G4bool SetMaximumEpsilonStep( G4double newEpsMax );
 
    /**
     * Methods to set/get flag for field changing energy.
     * For electric field this should be true; for magnetic field this
     * should be false.
     */
    inline G4bool DoesFieldChangeEnergy() const;
    inline void SetFieldChangesEnergy(G4bool value);
    
    /**
     * Needed for multi-threading, create and returns an allocated clone
     * of this object.
     */
    virtual G4FieldManager* Clone() const;

    /**
     * Static methods to set/get the maximum accepted epsilon.
     * If setting fails, with softFail=true it gives Warning, else
     * a FatalException.
     */
    static G4double GetMaxAcceptedEpsilon();
    static G4bool SetMaxAcceptedEpsilon(G4double maxEps, G4bool softFail= false);
   
  protected:

    /**
     * Logger for reporting on correctness of the proposed epsilon value.
     */
    void ReportBadEpsilonValue(G4ExceptionDescription& erm, G4double value,
                               const G4String& name) const;

    /** Epsilon_min/max values must be smaller than this for robust integration. */
    static G4double fMaxAcceptedEpsilon;
    static constexpr G4double fMinAcceptedEpsilon = 1000.0 * std::numeric_limits<G4double>::epsilon();

    /** Setting larger value will give warning. */
    static constexpr G4double fMaxWarningEpsilon = 0.001;

    /** Will not accept larger values. */
    static constexpr G4double fMaxFinalEpsilon = 0.02;
   
    /** Controls verbosity of constructors. */
    static G4bool fVerboseConstruction;
  
  private:

    /**
     * Checks whether the field/equation changes the energy and sets the data
     * member accordingly. Does not handle special cases - this must be done
     * separately (e.g. magnetic monopole in B field).
     */
    void InitialiseFieldChangesEnergy();
  
  private:

    /** Dependent objects -- with state that depends on tracking. */
    G4Field* fDetectorField = nullptr;
    G4ChordFinder* fChordFinder = nullptr;

    /** Flag to indicate if "new" was used to create the Chord Finder. */
    G4bool fAllocatedChordFinder = false; // 

    //  1. CHARACTERISTIC of field

    G4bool fFieldChangesEnergy = false;

    //  2. PARAMETERS that determine the accuracy of integration or intersection

    /** Value for the required accuracies for one tracking/physics step. */
    G4double fDelta_One_Step_Value = G4FieldDefaults::kDeltaOneStep;

    /** Value for the required accuracies for boundary intersection. */
    G4double fDelta_Intersection_Val = G4FieldDefaults::kDeltaIntersection;

    /** Values for the small possible relative accuracy of a step
        (corresponding to the greatest possible integration accuracy). */
    G4double fEpsilonMin = G4FieldDefaults::kMinimumEpsilonStep;
    G4double fEpsilonMax = G4FieldDefaults::kMaximumEpsilonStep;

    /** Global field manager set by G4TransportationManager to allow accessing
        the global field without dependency on navigation. */
    static G4ThreadLocal G4FieldManager* fGlobalFieldManager;
};

// Implementation of inline functions

#include "G4FieldManager.icc"

#endif

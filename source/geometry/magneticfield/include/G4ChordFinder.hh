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
// G4ChordFinder
//
// Class description:
//
// A class that provides RK integration of motion ODE (as does g4magtr)
// and also has a method that returns an Approximate point on the curve 
// near to a (chord) point.

// Author: John Apostolakis (CERN), 25.02.1997 - Design and implementation
// -------------------------------------------------------------------
#ifndef G4CHORDFINDER_HH
#define G4CHORDFINDER_HH

#include "G4VIntegrationDriver.hh"
#include "G4FieldParameters.hh"
#include "G4MagIntegratorStepper.hh"

#include <memory>

class G4VFSALIntegrationStepper;

class G4MagneticField;
class G4CachedMagneticField;
class G4HelixHeum;
class G4QSStepper;

/**
 * @brief G4ChordFinder is a class that provides Runge-Kutta integration of
 * motion ODE and also has a method that returns an approximate point on the
 * curve near to a (chord) point.
 */

class G4ChordFinder
{
  public:

    enum kIntegrationType {kDefaultDriverType=0, kFSALStepperType=1, 
                           kTemplatedStepperType, kRegularStepperType,
                           kBfieldDriverType, kQss2DriverType, kQss3DriverType};

    /**
     * The most flexible constructor, which allows the user to specify
     * any type of field, equation, stepper and integration driver.
     *  @param[in] pIntegrationDriver Pointer to the integrator driver to use.
     */
    explicit G4ChordFinder( G4VIntegrationDriver* pIntegrationDriver );

    /**
     * Constructor that creates defaults for all "children" classes.
     * The type of equation of motion is fixed.
     * A default type of stepper (Dormand Prince since release 10.4) is used,
     * and the corresponding integration driver.
     *  @param[in] itsMagField Pointer to the magnetic field.
     *  @param[in] stepMinimum Pointer to the magnetic field.
     *  @param[in] pItsStepper Optional pointer to the stepper algorithm.
     *  @param[in] stepperDriverChoice Type of stepper driver.
     */
    G4ChordFinder( G4MagneticField* itsMagField,
                   G4double stepMinimum = G4FieldDefaults::kMinimumStep,
                   G4MagIntegratorStepper* pItsStepper = nullptr,
                   G4int stepperDriverChoice = kTemplatedStepperType );

    /**
     * Destructor.
     */
    ~G4ChordFinder();

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4ChordFinder(const G4ChordFinder&) = delete;
    G4ChordFinder& operator=(const G4ChordFinder&) = delete;

    /**
     * Computes the step to take, based on chord limits.
     * Uses ODE solver's driver to find the endpoint that satisfies 
     * the chord criterion that: d_chord < delta_chord.
     *  @param[in,out] yCurrent The current track in field.
     *  @param[in] stepInitial Proposed initial step length.
     *  @param[in] epsStep_Relative Requested accuracy.
     *  @param[in] latestSafetyOrigin Last safety origin point. Unused.
     *  @param[in] lasestSafetyRadius Last safety distance. Unused.
     *  @returns The length of step taken.
     */
    inline G4double AdvanceChordLimited( G4FieldTrack& yCurrent,
                                         G4double stepInitial,
                                         G4double epsStep_Relative,
                                   const G4ThreeVector& latestSafetyOrigin,
                                         G4double lasestSafetyRadius );
     
    /**
     * Uses the Brent algorithm when possible, to determine the closest point
     * on the curve. Given a starting curve point A (CurveA_PointVelocity),
     * curve point B (CurveB_PointVelocity), a point E which is (generally)
     * not on the curve and  a point F which is on the curve (first
     * approximation), find new point S on the curve closer to point E. 
     * While advancing towards S utilise 'eps_step' as a measure of the
     * relative accuracy of each Step.
     *  @returns The end point on the curve closer to the given point E.
     */
    G4FieldTrack ApproxCurvePointS( const G4FieldTrack&  curveAPointVelocity,
                                    const G4FieldTrack&  curveBPointVelocity,
                                    const G4FieldTrack&  ApproxCurveV,
                                    const G4ThreeVector& currentEPoint,
                                    const G4ThreeVector& currentFPoint,
                                    const G4ThreeVector& PointG,
                                          G4bool first,  G4double epsStep );
 
    /**
     * If r=|AE|/|AB|, and s=true path lenght (AB)
     * returns the point that is r*s along the curve.
     */
    G4FieldTrack ApproxCurvePointV( const G4FieldTrack&  curveAPointVelocity,
                                    const G4FieldTrack&  curveBPointVelocity,
                                    const G4ThreeVector& currentEPoint,
                                          G4double       epsStep);

    /**
     * Calculates the inverse parabolic through the three points (x,y) and
     * returns the value x that, for the inverse parabolic, corresponds to y=0.
     */
    inline G4double InvParabolic( const G4double xa, const G4double ya,
                                  const G4double xb, const G4double yb,
                                  const G4double xc, const G4double yc );

    /**
     * Accessors and modifiers.
     */
    inline G4double GetDeltaChord() const;
    inline void SetDeltaChord(G4double newval);
    inline void SetIntegrationDriver(G4VIntegrationDriver* IntegrationDriver);
    inline G4VIntegrationDriver* GetIntegrationDriver();

    /**
     * Clears the internal state (last step estimate).
     */
    inline void ResetStepEstimate();

    /**
     * Sets the verbosity.
     *  @returns The old verbosity value.
     */
    inline G4int SetVerbose( G4int newvalue=1 ); 

    /**
     * Dispatch interface method for computing step.
     */
    inline void OnComputeStep(const G4FieldTrack* track);

    /**
     * Writes out to stream the parameters/state of the driver.
     */
    friend std::ostream& operator<<( std::ostream& os, const G4ChordFinder& cf);

    /**
     * Sets verbosity for constructor.
     */
    static void SetVerboseConstruction(G4bool v = true);

  private:  // ............................................................

    static G4bool gVerboseCtor;  // Verbosity for contructor

    //  Constants
    //  ---------------------
    const G4double fDefaultDeltaChord = G4FieldDefaults::kDeltaChord;

    //  PARAMETERS 
    //  ---------------------
    G4double  fDeltaChord;               //  Maximum miss distance 

    G4int fStatsVerbose = 0;  // if > 0, print Statistics in destructor

    //  DEPENDENT Objects
    //  ---------------------
    G4VIntegrationDriver*      fIntgrDriver = nullptr;
    G4MagIntegratorStepper*    fRegularStepperOwned = nullptr;
    G4MagIntegratorStepper*    fNewFSALStepperOwned = nullptr;
    std::unique_ptr<G4HelixHeum> fLongStepper;
    G4CachedMagneticField*     fCachedField = nullptr;
    G4QSStepper*               fQssStepperOwned = nullptr;
    G4EquationOfMotion*        fEquation = nullptr;  
};

// Inline function implementation:

#include "G4ChordFinder.icc"

#endif

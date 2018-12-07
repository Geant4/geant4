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
//
// 
// Class G4ChordFinder
//
// class description:
//
// A class that provides RK integration of motion ODE  (as does g4magtr)
// and also has a method that returns an Approximate point on the curve 
// near to a (chord) point.

// History:
// - 25.02.97 - John Apostolakis - Design and implementation 
// -------------------------------------------------------------------

#ifndef G4CHORDFINDER_HH
#define G4CHORDFINDER_HH

#include "G4VIntegrationDriver.hh"
#include "G4MagIntegratorStepper.hh"
class G4VFSALIntegrationStepper;
// #include "G4VFSALIntegratorSteper.hh"

class G4MagneticField;
class G4CachedMagneticField;

class G4ChordFinder
{
   public:  // with description

      explicit G4ChordFinder( G4VIntegrationDriver* pIntegrationDriver );
        // The most flexible constructor, which allows the user to specify
        //   any type of field, equation, stepper and integration driver.

      G4ChordFinder( G4MagneticField* itsMagField,
                     G4double         stepMinimum = 1.0e-2, // * mm 
                     G4MagIntegratorStepper* pItsStepper = nullptr,
                     // G4bool           useHigherEfficiencyStepper = true,
                     G4bool           useFSALstepper = false  );
        // A constructor that creates defaults for all "children" classes.
        //
        //   The type of equation of motion is fixed.
        //   A default type of stepper (Dormand Prince since release 10.4) is used,
        //   and the corresponding (old-style without themplates) integration driver.
        //   Except if 'useFSAL' is set (true), which provides a FSAL stepper
        //    and its corresponding specialised (templated) driver.
      
      virtual ~G4ChordFinder();

      inline G4double AdvanceChordLimited( G4FieldTrack& yCurrent,
                                           G4double stepInitial,
                                           G4double epsStep_Relative,
                                           const G4ThreeVector& latestSafetyOrigin,
                                           G4double lasestSafetyRadius);
        // Uses ODE solver's driver to find the endpoint that satisfies 
        // the chord criterion: that d_chord < delta_chord
        // -> Returns Length of Step taken.
     
      G4FieldTrack ApproxCurvePointS( const G4FieldTrack&  curveAPointVelocity,
                                      const G4FieldTrack&  curveBPointVelocity,
                                      const G4FieldTrack&  ApproxCurveV,
                                      const G4ThreeVector& currentEPoint,
                                      const G4ThreeVector& currentFPoint,
                                      const G4ThreeVector& PointG,
                                            G4bool first,  G4double epsStep);
 
      G4FieldTrack ApproxCurvePointV( const G4FieldTrack&  curveAPointVelocity,
                                      const G4FieldTrack&  curveBPointVelocity,
                                      const G4ThreeVector& currentEPoint,
                                            G4double       epsStep);

      inline G4double InvParabolic( const G4double xa, const G4double ya,
                                    const G4double xb, const G4double yb,
                                    const G4double xc, const G4double yc );

      inline G4double  GetDeltaChord() const;
      inline void      SetDeltaChord(G4double newval);

      inline void SetIntegrationDriver(G4VIntegrationDriver* IntegrationDriver);
      inline G4VIntegrationDriver* GetIntegrationDriver();
        // Access and set Driver.

      inline void ResetStepEstimate();
        // Clear internal state (last step estimate)

        // A report with the above -- and possibly other stats
      inline G4int SetVerbose( G4int newvalue=1); 
        // Set verbosity and return old value

      void OnComputeStep();

   protected:   // .........................................................
      void     PrintDchordTrial(G4int     noTrials, 
                                G4double  stepTrial, 
                                G4double  oldStepTrial, 
                                G4double  dChordStep);

   private:  // ............................................................

      G4ChordFinder(const G4ChordFinder&);
      G4ChordFinder& operator=(const G4ChordFinder&);
        // Private copy constructor and assignment operator.

   private:  // ............................................................
                                          // G4int    nOK, nBAD;

      // Constants
      const G4double fDefaultDeltaChord;  // SET in G4ChordFinder.cc = 0.25 mm

      //  PARAMETERS 
      //  ---------------------
      G4double  fDeltaChord;               //  Maximum miss distance 

      // G4double  fMultipleRadius;        //  Use forgotten
      G4int     fStatsVerbose;  // if > 0, print Statistics in destructor

      //  DEPENDENT Objects
      //  ---------------------
      G4VIntegrationDriver*      fIntgrDriver;
      G4MagIntegratorStepper*    fRegularStepperOwned= nullptr;
      G4MagIntegratorStepper*    fNewFSALStepperOwned= nullptr;
      G4CachedMagneticField*     fCachedField= nullptr;
   // G4VFSALIntegrationStepper* fOldFSALStepperOwned= nullptr;
      G4EquationOfMotion*        fEquation;  
};

// Inline function implementation:

#include "G4ChordFinder.icc"

#endif  // G4CHORDFINDER_HH

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
// $Id: G4ChordFinder.hh,v 1.21 2008-10-29 14:17:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include "G4MagIntegratorDriver.hh"
#include "G4FieldTrack.hh"

class G4MagneticField;  

class G4ChordFinder
{ 
   public:  // with description

      G4ChordFinder( G4MagInt_Driver* pIntegrationDriver );

      G4ChordFinder( G4MagneticField* itsMagField,
                     G4double         stepMinimum = 1.0e-2 * mm, 
                     G4MagIntegratorStepper* pItsStepper = 0 );  
        // A constructor that creates defaults for all "children" classes.
      
      virtual ~G4ChordFinder();

      G4double    AdvanceChordLimited( G4FieldTrack& yCurrent,
                                       G4double stepInitial,
                                       G4double epsStep_Relative,
                                       const G4ThreeVector latestSafetyOrigin,
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

      inline void SetChargeMomentumMass(G4double pCharge,  // in e+ units
                                        G4double pMomentum,
                                        G4double pMass );
        // Function to inform integration driver of charge, speed.

      inline void SetIntegrationDriver(G4MagInt_Driver* IntegrationDriver);
      inline G4MagInt_Driver* GetIntegrationDriver();
        // Access and set Driver.

      inline void ResetStepEstimate();
        // Clear internal state (last step estimate)

      inline G4int GetNoCalls(); 
      inline G4int GetNoTrials();        // Total number of trials
      inline G4int GetNoMaxTrials();     // Maximum # of trials for one call
        // Get statistics about number of calls & trials in FindNextChord

      virtual void   PrintStatistics(); 
        // A report with the above -- and possibly other stats
      inline G4int SetVerbose( G4int newvalue=1); 
        // Set verbosity and return old value

      void SetFractions_Last_Next( G4double fractLast= 0.90, 
                                   G4double fractNext= 0.95 ); 
        // Parameters for  performance ... change with great care

      inline void SetFirstFraction(G4double fractFirst);
        // Parameter for  performance ... change with great care

   public:  // without description

      void     TestChordPrint( G4int    noTrials, 
                               G4int    lastStepTrial, 
                               G4double dChordStep, 
                               G4double nextStepTrial );

        //   Printing for monitoring ...
 
      inline   G4double GetFirstFraction();         // Originally 0.999
      inline   G4double GetFractionLast();          // Originally 1.000
      inline   G4double GetFractionNextEstimate();  // Originally 0.980
      inline   G4double GetMultipleRadius();        // No original value
        //  Parameters for adapting performance ... use with great care

   protected:   // .........................................................

      inline  void    AccumulateStatistics( G4int noTrials ); 
        // Accumulate the basic statistics 
        //   - other specialised ones must be kept by derived classes
 
      inline G4bool AcceptableMissDist(G4double dChordStep) const;

      G4double NewStep( G4double stepTrialOld, 
                        G4double dChordStep,     // Current dchord estimate
                        G4double& stepEstimate_Unconstrained ) ;  
      
      virtual G4double FindNextChord( const  G4FieldTrack& yStart,
                              G4double     stepMax,
                              G4FieldTrack& yEnd,
                              G4double&    dyErr,      //  Error of endpoint 
                              G4double     epsStep,
                              G4double*  pNextStepForAccuracy,  // = 0,
                              const G4ThreeVector latestSafetyOrigin,
                              G4double       latestSafetyRadius 
                                      );  

      void     PrintDchordTrial(G4int     noTrials, 
                                G4double  stepTrial, 
                                G4double  oldStepTrial, 
                                G4double  dChordStep);

      inline G4double GetLastStepEstimateUnc(); 
      inline void     SetLastStepEstimateUnc( G4double stepEst ); 

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
      //    Internal parameters
      G4double  fFirstFraction, fFractionLast, fFractionNextEstimate;
      G4double  fMultipleRadius; 
      G4int     fStatsVerbose;  // if > 0, print Statistics in destructor

      //  DEPENDENT Objects
      //  ---------------------
      G4MagInt_Driver*        fIntgrDriver;
      G4MagIntegratorStepper* fDriversStepper; 
      G4bool                  fAllocatedStepper;  // Bookkeeping of dependent object
      G4EquationOfMotion*     fEquation; 

      //  STATE information
      //  --------------------
      G4double    fLastStepEstimate_Unconstrained;
        //  State information for efficiency

      // For Statistics
      // -- G4int   fNoTrials, fNoCalls;
      G4int   fTotalNoTrials_FNC,  fNoCalls_FNC, fmaxTrials_FNC; // fnoTimesMaxTrFNC; 
};

// Inline function implementation:

#include "G4ChordFinder.icc"

#endif  // G4CHORDFINDER_HH

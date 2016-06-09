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
// $Id: G4ChordFinder.hh,v 1.16 2003/11/13 17:53:47 japost Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// class G4ChordFinder
//
// Class description:
//
// A class that provides RK integration of motion ODE  (as does g4magtr)
// and also has a method that returns an Approximate point on the curve 
// near to a (chord) point.

// History:
// - 25.02.97 John Apostolakis,  design and implementation 
// - 05.03.97 V. Grichine , makeup to G4 'standard'
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

      G4FieldTrack ApproxCurvePointV(const  G4FieldTrack&  curveAPointVelocity,
                                     const  G4FieldTrack&  curveBPointVelocity,
                                     const  G4ThreeVector& currentEPoint,
                                            G4double      epsStep);

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

   protected:   // .........................................................

      inline  void    AccumulateStatistics( G4int noTrials ); 
        // Accumulate the basic statistics 
        //   - other specialised ones must be kept by derived classes
 
      inline G4bool AcceptableMissDist(G4double dChordStep) const;

      G4double NewStep( G4double stepTrialOld, 
                        G4double dChordStep,     // Current dchord estimate
                        G4double& stepEstimate_Unconstrained ) ;  
      
      virtual G4double FindNextChord( const  G4FieldTrack  yStart,
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
  public:  // no description 
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

   public:  // with description 
      void     SetFractions_Last_Next( G4double fractLast= 0.90, 
				       G4double fractNext= 0.95 ); 
      //  Parameters for  performance ... change with great care

      inline   void     SetFirstFraction(G4double fractFirst);
      //  Parameter for  performance ... change with great care

   protected:
      inline G4double GetLastStepEstimateUnc(); 
      inline void     SetLastStepEstimateUnc( G4double stepEst ); 

   private:  // ............................................................

      G4ChordFinder(const G4ChordFinder&);
      G4ChordFinder& operator=(const G4ChordFinder&);
        // Private copy constructor and assignment operator.

   private:  // ............................................................
                                            // G4int    nOK, nBAD;
      G4MagInt_Driver* fIntgrDriver;

      const G4double fDefaultDeltaChord;  // SET in G4ChordFinder.cc = 0.25 mm

      G4double fDeltaChord;                        //  Maximum miss distance 

      G4double    fLastStepEstimate_Unconstrained; //  State information for efficiency
      //  Variables used in construction/destruction
      G4bool fAllocatedStepper;
      G4EquationOfMotion* fEquation; 
      G4MagIntegratorStepper* fDriversStepper; 

      //  Parameters 
      G4double  fFirstFraction, fFractionLast, fFractionNextEstimate;
      G4double  fMultipleRadius; 

      // For Statistics
      // -- G4int   fNoTrials, fNoCalls;
      G4int   fTotalNoTrials_FNC,  fNoCalls_FNC, fmaxTrials_FNC; // fnoTimesMaxTrFNC; 
      G4int   fStatsVerbose;  // if > 0, print Statistics in destructor
};

// Inline function implementation:

#include "G4ChordFinder.icc"

#endif  // G4CHORDFINDER_HH

// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ChordFinder.hh,v 1.6 2001-03-23 18:50:33 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#ifndef G4CHORDFINDER_HH
#define G4CHORDFINDER_HH

#include "G4MagIntegratorDriver.hh"
#include "G4FieldTrack.hh"
#include "G4MagneticField.hh"
                                             // class G4Mag_EqRhs;  
                                             // class G4MagIntegratorStepper; 

class G4ChordFinder
{ 
   public:  // with description

      G4ChordFinder( G4MagInt_Driver* pIntegrationDriver );

      G4ChordFinder( G4MagneticField* itsMagField,
		     G4double         stepMinimum = 1.0e-2 * mm, 
		     G4MagIntegratorStepper* pItsStepper = 0 );  
        // A constructor that creates defaults for all "children" classes.
                  
      ~G4ChordFinder();


      G4double    AdvanceChordLimited( G4FieldTrack& yCurrent,
                                       G4double stepInitial,
                                       G4double epsStep_Relative  );
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

   protected:   // .........................................................

      inline G4bool AcceptableMissDist(G4double dChordStep) const;

      G4double NewStep( G4double stepTrialOld, 
		        G4double dChordStep,     // Current dchord estimate
                        G4double& stepEstimate_Unconstrained ) ;  
      
      G4double FindNextChord( const  G4FieldTrack  yStart,
		              G4double     stepMax,
			      G4FieldTrack& yEnd,
			      G4double&    dyErr,      //  Error of endpoint 
			      G4double     epsStep );  

   private:  // ............................................................

      G4ChordFinder(const G4ChordFinder&);
      G4ChordFinder& operator=(const G4ChordFinder&);
        // Private copy constructor and assignment operator.

   private:  // ............................................................
                                            // G4int    nOK, nBAD;
      G4MagInt_Driver* fIntgrDriver;

      G4double fDeltaChord;                        //  Maximum miss distance 

      G4double    fLastStepEstimate_Unconstrained; //  State information for efficiency

      static const G4double fDefaultDeltaChord;  // SET in G4ChordFinder.cc = 3 mm

      //  Variables used in construction/destruction
      G4bool fAllocatedStepper;
      G4EquationOfMotion* fEquation; 
      G4MagIntegratorStepper* fDriversStepper; 
};

// Inline function implementation:

#include "G4ChordFinder.icc"

#endif  // G4CHORDFINDER_HH

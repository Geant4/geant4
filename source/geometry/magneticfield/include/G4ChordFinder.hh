// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ChordFinder.hh,v 1.1 1999-01-07 16:07:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------------------
//	GEANT 4  include file implementation
//
//	For information related to this code contact:
//	CERN, IT Division (formely CN), ASD group
// ------------------------------------------------------------------------
//
//  A class that provides RK integration of motion ODE  (as does g4magtr)
//   and also has a method that returns an Approximate point on the curve 
//   near to a (chord) point.
//
// 25.02.97 John Apostolakis,  design and implementation 
// 05.03.97 V. Grichine , makeup to G4 'standard'

#ifndef G4CHORDFINDER_HH
#define G4CHORDFINDER_HH

                                             // #include "globals.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4FieldTrack.hh"
#include "G4MagneticField.hh"
                                             // class G4Mag_EqRhs;  
                                             // class G4MagIntegratorStepper; 

class G4ChordFinder
{ 
   public:                      // Constructors 

      G4ChordFinder( G4MagInt_Driver* pIntegrationDriver );

      // A constructor that creates defaults for all "children" classes
      //
      G4ChordFinder( G4MagneticField* itsMagField,
		     G4double         stepMinimum = 1.0e-2 * mm, 
		     G4MagIntegratorStepper* pItsStepper = 0 );  
                  
      ~G4ChordFinder();

      //      Uses ODE solver's driver to find the endpoint that satisfies 
      //   the chord criterion: that d_chord < delta_chord
      //   -> Returns Length of Step taken

      G4double    AdvanceChordLimited( G4FieldTrack& yCurrent,
                                    const  G4double     stepInitial,
                                    const  G4double     epsStep  );

      G4FieldTrack ApproxCurvePointV(const  G4FieldTrack&  curveAPointVelocity,
		  	             const  G4FieldTrack&  curveBPointVelocity,
			             const  G4ThreeVector& currentEPoint,
			             const  G4double      epsStep);

      G4double  GetDeltaChord();
      void      SetDeltaChord( G4double newval);

      // Routine to inform integration driver of charge, speed 
      //  
      void SetChargeMomentumMass( const G4double pCharge,    // in e+ units
				  const G4double pMomentum,
				  const G4double pMass );

      // Access and set Driver
      //
      void             SetIntegrationDriver( G4MagInt_Driver* IntegrationDriver)
                                          { fIntgrDriver=IntegrationDriver;}
      G4MagInt_Driver* GetIntegrationDriver()
                                          { return fIntgrDriver;}
      
   protected:   // .........................................................

      G4bool AcceptableMissDist(G4double dChordStep) 
      { 
 	return (dChordStep <= fDeltaChord) ;
      }

      G4double NewStep( const G4double stepTrialOld, 
		        const G4double dChordStep ) ;  // Current dchord 
      
      G4double FindNextChord( const  G4FieldTrack  yStart,
		              const  G4double     stepMax,
			      G4FieldTrack& yEnd,
			      G4double&    dyErr,      //  Error of endpoint 
			      G4double     epsStep );  

   private:  // ............................................................
                                            // G4int    nOK, nBAD;
      G4MagInt_Driver* fIntgrDriver;

      G4double fDeltaChord;    

      static const G4double fDefaultDeltaChord;  // SET in G4ChordFinder.cc = 3 mm

      //  Variables used in construction/destruction
      G4bool fAllocatedStepper;
      G4Mag_EqRhs* fEquation; 
      G4MagIntegratorStepper* fDriversStepper; 
};



// Inline function implementation:

inline
G4ChordFinder:: G4ChordFinder( G4MagInt_Driver* pIntegrationDriver )
: fDeltaChord( fDefaultDeltaChord )
{
    fIntgrDriver= pIntegrationDriver ;
    fAllocatedStepper= false ;
} 

inline void
G4ChordFinder::SetChargeMomentumMass( const G4double pCharge,  // in e+ units
				      const G4double pMomentum,
				      const G4double pMass )
{
   fIntgrDriver-> SetChargeMomentumMass(pCharge, pMomentum, pMass);
}

inline G4double  G4ChordFinder::GetDeltaChord() 
{   return fDeltaChord; }

inline void      G4ChordFinder::SetDeltaChord( G4double newval)
{   fDeltaChord=newval; }

#endif  // G4CHORDFINDER_HH

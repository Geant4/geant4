// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FieldTrack.hh,v 1.2 1999-12-15 14:49:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//  Data structure bringing together a magnetic track's state.
//   (position, momentum direction & modulus, energy, spin, ... )
//  Uses/abilities:
//   - does not maintain any relationship between its data (eg energy/momentum)
//   - for use in Runge-Kutta solver (in passing it the values right now).
//
//  First version: Oct 14, 1996  John Apostolakis
//  Modified:      Oct 24, 1996  JA: Added dist_on_curve, deleted constructor
//                 Nov  5, 1998  JA: Added energy, momentum, TOF, spin
//                                   & several constructor, access, set methods
//                                   
// 
#ifndef G4FieldTrack_HH
#define G4FieldTrack_HH

#include "G4ThreeVector.hh"

class  G4FieldTrack{
   public:
     //  Constructors

     G4FieldTrack( const G4ThreeVector& pPosition, 
		   const G4ThreeVector& pVelocity,   // Or UnitVelocity
		   const G4double       curve_length,
		   const G4double       Energy,
		   const G4double       LabratTimeOfFlight=0.0,
		   const G4double       ProperTimeOfFlight=0.0, 
		   const G4ThreeVector* pSpin=0);

     G4FieldTrack( const G4FieldTrack&   pFieldTrack );

     // Destructor 
     ~G4FieldTrack();

     // Equality operator
     G4FieldTrack& operator = ( const G4FieldTrack & rStVec );

     // Old multi-set method
     inline G4FieldTrack&  SetCurvePnt( 
				const G4ThreeVector& pPosition, 
				const G4ThreeVector& pVelocity,
				const G4double       s_curve );

     //  Access Methods:   ("Const")
     G4ThreeVector  Position() const;   // Renamed to GetPosition
     G4ThreeVector  GetVelocity() const;   
     G4double       CurveS()    const;  // distance along curve of point
     //   Old methods above to be deleted. 
     G4ThreeVector  GetPosition() const; 
     const G4ThreeVector& GetMomentumDir() const;
     G4double       GetCurveLength() const;  // distance along curve of point
     // G4double       GetEnergy() const;       //  Wrong Energy  --> FIXME
     G4double       GetMomentumModulus() const;
     G4ThreeVector  GetSpin()   const;
     G4double       GetLabTimeOfFlight() const;
     G4double       GetProperTimeOfFlight() const;

     //  Modifiers
     void SetPosition(G4ThreeVector nPos); 
     void SetVelocity(G4ThreeVector nMomDir);    // does change mom-dir too
     void SetMomentumDir(G4ThreeVector nMomDir); // does NOT change velocity
     void SetCurveLength(G4double nCurve_s); // distance along curve 
     void SetEnergy(G4double nEnergy);              // does not modify momentum
     void SetMomentumModulus(G4double nMomentumMod); // does not modify energy
     void SetSpin(G4ThreeVector nSpin);
     void SetLabTimeOfFlight(G4double nTOF); 
     void SetProperTimeOfFlight(G4double nTOF);
     // older one: 
     void SetCurveS(G4double new_curve_s);

     // G4double*      PosVelVec();       // [6]  Needed for RK integrator
     //  This old method completely broke encapsulation ?  
  
     //  Needed and should be used only for RK integration driver
     // static const G4int ncompSVEC=15;
     enum { ncompSVEC = 16 };
     void DumpToArray(   G4double valArr[ncompSVEC] ) const; 
     void LoadFromArray( const G4double valArr[ncompSVEC] ); 
     
     friend  G4std::ostream& operator<<( G4std::ostream& os, G4FieldTrack& SixVec);

   private:
     G4double  SixVector[6];
     G4double  fDistanceAlongCurve;  // distance along curve of point
     G4double  fEnergy;
     G4double  fLabTimeOfFlight;
     G4double  fProperTimeOfFlight;
     G4double  fMomentumModulus;
     G4ThreeVector fSpin;
     G4ThreeVector fMomentumDir;
}; 

#include "G4FieldTrack.icc"

#endif  /* End of ifndef G4FieldTrack_HH */

// Rename:
//
// s/distance_along_curve/fDistanceAlongCurve/g;
// s/SixVector/fSixVector/g;
// s/G4SixVector/G4FieldTrack/g;

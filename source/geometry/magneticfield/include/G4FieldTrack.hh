// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FieldTrack.hh,v 1.3 2000-04-27 09:14:05 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4FieldTrack
//
// Class description:
//
// Data structure bringing together a magnetic track's state.
// (position, momentum direction & modulus, energy, spin, ... )
// Uses/abilities:
//  - does not maintain any relationship between its data (eg energy/momentum).
//  - for use in Runge-Kutta solver (in passing it the values right now).

// History
// - First version: Oct 14, 1996  John Apostolakis
// - Modified:      Oct 24, 1996  JA: Added dist_on_curve, deleted constructor
//                  Nov  5, 1998  JA: Added energy, momentum, TOF, spin &
//                                    several constructor, access, set methods

#ifndef G4FieldTrack_HH
#define G4FieldTrack_HH

#include "G4ThreeVector.hh"

class  G4FieldTrack
{
   public:  // with description

     G4FieldTrack( const G4ThreeVector& pPosition, 
		   const G4ThreeVector& pVelocity,   // Or UnitVelocity
		   const G4double       curve_length,
		   const G4double       Energy,
		   const G4double       LabratTimeOfFlight=0.0,
		   const G4double       ProperTimeOfFlight=0.0, 
		   const G4ThreeVector* pSpin=0);

     G4FieldTrack( const G4FieldTrack&   pFieldTrack );

     ~G4FieldTrack();
       // Destructor 

     G4FieldTrack& operator = ( const G4FieldTrack & rStVec );
       // Equality operator

     G4ThreeVector  GetVelocity() const;   
     G4ThreeVector  GetPosition() const; 
     const G4ThreeVector& GetMomentumDir() const;
     G4double       GetCurveLength() const;
       // Distance along curve of point.
     G4double       GetMomentumModulus() const;
     G4ThreeVector  GetSpin()   const;
     G4double       GetLabTimeOfFlight() const;
     G4double       GetProperTimeOfFlight() const;
       // Accessors.

     void SetPosition(G4ThreeVector nPos); 
     void SetVelocity(G4ThreeVector nMomDir);
       // Does change mom-dir too.
     void SetMomentumDir(G4ThreeVector nMomDir);
       // Does NOT change velocity.
     void SetCurveLength(G4double nCurve_s);
       // Distance along curve.
     void SetEnergy(G4double nEnergy);
       // Does not modify momentum.
     void SetMomentumModulus(G4double nMomentumMod);
       // Does not modify energy.
     void SetSpin(G4ThreeVector nSpin);
     void SetLabTimeOfFlight(G4double nTOF); 
     void SetProperTimeOfFlight(G4double nTOF);
       //  Modifiers

   public: // without description

     G4FieldTrack&  SetCurvePnt(const G4ThreeVector& pPosition, 
				const G4ThreeVector& pVelocity,
				const G4double       s_curve );
       // Old multi-set method

     G4ThreeVector  Position() const;
       // Renamed to GetPosition

     G4double       CurveS()    const;  // distance along curve of point
     void SetCurveS(G4double new_curve_s);
       // Old methods to be deleted. 

     // G4double       GetEnergy() const;       //  Wrong Energy  --> FIXME

     // G4double*      PosVelVec();       // [6]  Needed for RK integrator
       // This old method completely broke encapsulation ?  
  
     // static const G4int ncompSVEC=15;
       // Needed and should be used only for RK integration driver
     enum { ncompSVEC = 16 };
     void DumpToArray(   G4double valArr[ncompSVEC] ) const; 
     void LoadFromArray( const G4double valArr[ncompSVEC] ); 
     
     friend  G4std::ostream&
             operator<<( G4std::ostream& os, G4FieldTrack& SixVec);

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

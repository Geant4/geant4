// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FieldTrack.hh,v 1.4 2000-11-01 15:15:48 gcosmo Exp $
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
		         G4double       curve_length,
		         G4double       Energy,
		         G4double       LabratTimeOfFlight=0.0,
		         G4double       ProperTimeOfFlight=0.0, 
		   const G4ThreeVector* pSpin=0);

     G4FieldTrack( const G4FieldTrack&   pFieldTrack );

     ~G4FieldTrack();
       // Destructor 

     inline G4FieldTrack& operator = ( const G4FieldTrack & rStVec );
       // Assignment operator

     inline G4ThreeVector  GetVelocity() const;   
     inline G4ThreeVector  GetPosition() const; 
     inline const G4ThreeVector& GetMomentumDir() const;
     inline G4double       GetCurveLength() const;
       // Distance along curve of point.
     inline G4double       GetMomentumModulus() const;
     inline G4ThreeVector  GetSpin()   const;
     inline G4double       GetLabTimeOfFlight() const;
     inline G4double       GetProperTimeOfFlight() const;
       // Accessors.

     inline void SetPosition(G4ThreeVector nPos); 
     inline void SetVelocity(G4ThreeVector nMomDir);
       // Does change mom-dir too.
     inline void SetMomentumDir(G4ThreeVector nMomDir);
       // Does NOT change velocity.
     inline void SetCurveLength(G4double nCurve_s);
       // Distance along curve.
     inline void SetEnergy(G4double nEnergy);
       // Does not modify momentum.
     inline void SetMomentumModulus(G4double nMomentumMod);
       // Does not modify energy.
     inline void SetSpin(G4ThreeVector nSpin);
     inline void SetLabTimeOfFlight(G4double nTOF); 
     inline void SetProperTimeOfFlight(G4double nTOF);
       //  Modifiers

   public: // without description

     inline G4FieldTrack& SetCurvePnt(const G4ThreeVector& pPosition, 
				      const G4ThreeVector& pVelocity,
				            G4double       s_curve );
       // Old multi-set method

     inline G4ThreeVector Position() const;
       // Renamed to GetPosition

     inline G4double CurveS() const;
       // Distance along curve of point
     inline void SetCurveS(G4double new_curve_s);
       // Old methods to be deleted. 

     // G4double       GetEnergy() const;       //  Wrong Energy  --> FIXME

     // G4double*      PosVelVec();       // [6]  Needed for RK integrator
       // This old method completely broke encapsulation ?  
  
     // static const G4int ncompSVEC=15;
       // Needed and should be used only for RK integration driver
     enum { ncompSVEC = 16 };
     inline void DumpToArray(G4double valArr[ncompSVEC]) const; 
     inline void LoadFromArray(G4double valArr[ncompSVEC]); 
     
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

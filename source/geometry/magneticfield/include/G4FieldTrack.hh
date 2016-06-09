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
// $Id: G4FieldTrack.hh,v 1.21 2006-11-13 18:24:35 gcosmo Exp $
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
//                  May 10, 2006  JA: Added charge, "default" constructor
// -------------------------------------------------------------------

#ifndef G4FieldTrack_HH
#define G4FieldTrack_HH

#include "G4ThreeVector.hh"

class  G4FieldTrack
{
   public:  // with description

     G4FieldTrack( const G4ThreeVector& pPosition, 
                         G4double       LaboratoryTimeOfFlight,
                   const G4ThreeVector& pMomentumDirection,
                         G4double       kineticEnergy,
                         G4double       restMass_c2,
                         G4double       charge, 
                   const G4ThreeVector& pSpin,
                         G4double       magnetic_dipole_moment= 0.0,
                         G4double       curve_length= 0.0
                 );
     G4FieldTrack( const G4FieldTrack&   pFieldTrack ); 
     G4FieldTrack( char );   //  Almost default constructor

     ~G4FieldTrack();
       // End of preferred Constructors / Destructor 

     inline void
     UpdateState( const G4ThreeVector& pPosition, 
                        G4double       LaboratoryTimeOfFlight,
                  const G4ThreeVector& pMomentumDirection,
                        G4double       kineticEnergy); 
        //  Update four-vectors for space/time and momentum/energy
        //    Also resets curve length.
     inline
     void  UpdateFourMomentum( G4double             kineticEnergy, 
                               const G4ThreeVector& momentumDirection ); 
        //  Update momentum (and direction), and kinetic energy 

     void SetChargeAndMoments(G4double charge, 
                              G4double magnetic_dipole_moment= DBL_MAX,
                              G4double electric_dipole_moment= DBL_MAX,
                              G4double magnetic_charge=DBL_MAX );
        //   Sets all charges and moments

     G4FieldTrack( const G4ThreeVector& pPosition, 
                   const G4ThreeVector& pMomentumDirection,
                         G4double       curve_length,
                         G4double       kineticEnergy,
                   const G4double       restMass_c2,
                         G4double       velocity,
                         G4double       LaboratoryTimeOfFlight=0.0,
                         G4double       ProperTimeOfFlight=0.0, 
                   const G4ThreeVector* pSpin=0);
          // Older constructor
          //  --->  Misses charge !!!

     inline G4FieldTrack& operator = ( const G4FieldTrack & rStVec );
       // Assignment operator

     inline G4ThreeVector  GetMomentum() const;   
     inline G4ThreeVector  GetPosition() const; 
     inline const G4ThreeVector& GetMomentumDir() const;
     inline G4ThreeVector  GetMomentumDirection() const;
     inline G4double       GetCurveLength() const;
       // Distance along curve of point.

     inline G4ThreeVector  GetSpin()   const;
     inline G4double       GetLabTimeOfFlight() const;
     inline G4double       GetProperTimeOfFlight() const;
     inline G4double       GetKineticEnergy() const;
     inline G4double       GetCharge() const;
       // Accessors.

     inline void SetPosition(G4ThreeVector nPos); 
     inline void SetMomentum(G4ThreeVector nMomDir);
       // Does change mom-dir too.

     inline void SetMomentumDir(G4ThreeVector nMomDir);
       // Does NOT change Momentum or Velocity Vector.

     inline void SetCurveLength(G4double nCurve_s);
       // Distance along curve.
     inline void SetKineticEnergy(G4double nEnergy);
       // Does not modify momentum.

     inline void SetRestMass(G4double Mass_c2) { fRestMass_c2= Mass_c2; }
     inline void SetSpin(G4ThreeVector nSpin);
     inline void SetLabTimeOfFlight(G4double nTOF); 
     inline void SetProperTimeOfFlight(G4double nTOF);
       //  Modifiers

   public: // without description

     inline G4FieldTrack& SetCurvePnt(const G4ThreeVector& pPosition, 
                                      const G4ThreeVector& pMomentum,
                                            G4double       s_curve );
     inline void          InitialiseSpin( const G4ThreeVector& Spin );
       //  Used to update / initialise the state

     enum { ncompSVEC = 12 };
       // Needed and should be used only for RK integration driver

     inline void DumpToArray(G4double valArr[ncompSVEC]) const; 
     inline void LoadFromArray(const G4double valArr[ncompSVEC],
                                     G4int noVarsIntegrated);
     friend  std::ostream&
             operator<<( std::ostream& os, const G4FieldTrack& SixVec);

   private:

     G4double  SixVector[6];
     G4double  fDistanceAlongCurve;  // distance along curve of point
     G4double  fKineticEnergy;
     G4double  fRestMass_c2;
     G4double  fLabTimeOfFlight;
     G4double  fProperTimeOfFlight;
     G4ThreeVector fSpin;
     G4ThreeVector fMomentumDir;

   private:   //  Implementation detail -- daughter class

       class G4ChargeState  // Charge & moments
       {
         public:  // without description

           inline G4ChargeState(G4double charge,                       
                                G4double magnetic_dipole_moment= 0.0,  
                                G4double electric_dipole_moment= 0.0,  
                                G4double magnetic_charge= 0.0);  
           inline G4ChargeState( const G4ChargeState& right ); 

           inline void SetCharge(G4double charge){ fCharge= charge; }

           // Revise the charge (in units of the positron charge)
           // do not change moments

           void SetChargeAndMoments(G4double charge, 
                                    G4double magnetic_dipole_moment= DBL_MAX,
                                    G4double electric_dipole_moment= DBL_MAX,
                                    G4double magnetic_charge=DBL_MAX );
           //  Revise the charge and all moments

           G4double GetCharge() const { return fCharge; }  
           G4double GetMagneticDipoleMoment() const { return fMagn_dipole; } 
           G4double ElectricDipoleMoment() const { return fElec_dipole; } 
           G4double MagneticCharge() const { return fMagneticCharge; } 
 
         private:

           G4double fCharge; 
           G4double fMagn_dipole;
           G4double fElec_dipole;
           G4double fMagneticCharge;  // for magnetic monopole
       };

       G4ChargeState fChargeState;

    public: // Access

       const G4ChargeState* GetChargeState() const { return &fChargeState; } 
}; 

#include "G4FieldTrack.icc"

#endif  /* End of ifndef G4FieldTrack_HH */

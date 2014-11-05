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
// $Id: G4ChargeState.hh $
//
//
// class G4ChargeState
//
// Class description:
//
// History
// - First version: Apr 10, 2013  John Apostolakis, Peter Gumplinger
// - Modified:
// -------------------------------------------------------------------

#ifndef G4ChargeState_HH
#define G4ChargeState_HH

#include "globals.hh"

class  G4ChargeState  // Charge & moments
{
   public:  // without description

     inline G4ChargeState(G4double charge,
                          G4double magnetic_dipole_moment,
                          G4double pdgSpin, 
                          G4double electric_dipole_moment = 0.0,
                          G4double magnetic_charge = 0.0);

     inline G4ChargeState( const G4ChargeState& right );
     inline G4ChargeState& operator = ( const G4ChargeState& right );

     void SetChargeSpinMoments(G4double charge,
                               G4double pdgSpin,  
                               G4double magnetic_dipole_moment= DBL_MAX,
                               G4double electric_dipole_moment= DBL_MAX,
                               G4double magnetic_charge= DBL_MAX );
     //  Revise the charge, pdgSpin, and optionally both moments and magnetic charge 

     void SetCharge(G4double charge){ fCharge = charge; }
     G4double GetCharge() const { return fCharge; }
     // Revise the charge (in units of the positron charge)


     //  Basic Get / Set methods 
     void     SetPDGSpin(G4double spin){ fSpin = spin; }
     G4double GetPDGSpin() const { return fSpin; }

     void     SetMagneticDipoleMoment(G4double moment){ fMagn_dipole = moment; }
     G4double GetMagneticDipoleMoment() const { return fMagn_dipole; }

     void     SetElectricDipoleMoment(G4double moment){ fElec_dipole = moment; }
     G4double ElectricDipoleMoment() const { return fElec_dipole; }

     void     SetMagneticCharge(G4double charge){ fMagneticCharge=charge; }
     G4double MagneticCharge() const { return fMagneticCharge; }

     // Auxiliary methods to set several properties at once 

     inline void SetChargeMdm(G4double charge, G4double mag_dipole_moment);
       // SetCharge and Magnetic Dipole Moment

     inline void SetChargeMdmSpin(G4double charge,
                                  G4double magnetic_dipole_moment,
                                  G4double pdgSpin); 

     inline void SetChargeSpin(G4double charge,
                               G4double pdgSpin); 

     //  Revise the charge, spin and all both moments
     inline void SetChargeDipoleMoments(G4double charge,
                                 G4double magnetic_dipole_moment,
                                 G4double electric_dipole_moment); 

     inline void SetChargesAndMoments(G4double charge,
                               G4double magnetic_dipole_moment, 
                               G4double electric_dipole_moment,
                               G4double magnetic_charge );
   



   public: // Obsolete
     inline void     SetSpin(G4double spin){ SetPDGSpin( spin); } 
     inline G4double GetSpin() const { return GetPDGSpin(); }

   private:

     G4double fCharge;
     G4double fSpin;
     G4double fMagn_dipole;
     G4double fElec_dipole;
     G4double fMagneticCharge;  // for magnetic monopole
};

// Inline methods implementation

inline G4ChargeState::G4ChargeState(G4double charge,
                             G4double magnetic_dipole_moment,
                             G4double spin,
                             G4double electric_dipole_moment,
                             G4double magnetic_charge)
{
   fCharge         = charge;
   fSpin           = spin;
   fMagn_dipole    = magnetic_dipole_moment;
   fElec_dipole    = electric_dipole_moment;
   fMagneticCharge = magnetic_charge;
}

inline G4ChargeState::G4ChargeState( const G4ChargeState& right )
{
  fCharge         = right.fCharge;
  fSpin           = right.fSpin;
  fMagn_dipole    = right.fMagn_dipole;
  fElec_dipole    = right.fElec_dipole;
  fMagneticCharge = right.fMagneticCharge;
}

inline G4ChargeState& G4ChargeState::operator = ( const G4ChargeState& right )
{
  if (&right == this) return *this;

  fCharge         = right.fCharge;
  fSpin           = right.fSpin;
  fMagn_dipole    = right.fMagn_dipole;
  fElec_dipole    = right.fElec_dipole;
  fMagneticCharge = right.fMagneticCharge;

  return *this;
}

inline void G4ChargeState::SetChargeMdm(G4double charge, G4double mag_dipole_moment)
{ 
   SetCharge( charge ); 
   SetMagneticDipoleMoment( mag_dipole_moment); 
} 

inline void G4ChargeState::SetChargeMdmSpin(G4double charge,
                                  G4double magDipoleMoment,
                                  G4double pdgSpin)
{
   SetChargeMdm( charge, magDipoleMoment ); 
   SetPDGSpin( pdgSpin ); 
}

inline void G4ChargeState::SetChargeSpin(G4double charge,
                                         G4double pdgSpin)
{
   SetCharge( charge ); 
   SetPDGSpin( pdgSpin ); 
}

inline void 
G4ChargeState::SetChargeDipoleMoments(G4double charge,
                                      G4double magneticDM,
                                      G4double electricDM)
{
   SetChargeMdm( charge, magneticDM ); 
   SetElectricDipoleMoment( electricDM ); 
}

inline void 
G4ChargeState::SetChargesAndMoments(G4double charge,
                                    G4double magneticDM,
                                    G4double electricDM,
                                    G4double magnetic_charge )
{
   SetChargeDipoleMoments( charge, magneticDM, electricDM); 
   SetMagneticCharge( magnetic_charge ); 
}
#endif  /* End of ifndef G4ChargeState_HH */

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
// G4ChargeState
//
// Class description:
//
// Container for magnetic charge and moments.

// Authors: J.Apostolakis (CERN), P.Gumplinger (TRIUMF), 10.04.2013  
// -------------------------------------------------------------------
#ifndef G4CHARGESTATE_HH
#define G4CHARGESTATE_HH

#include "globals.hh"

/**
 * @brief G4ChargeState is a container for magnetic charge and moments.
 */

class G4ChargeState
{
  public:

    /**
     * Constructor for G4ChargeState.
     *  @param[in] charge Particle charge.
     *  @param[in] magnetic_dipole_moment Magnetic dipole moment.
     *  @param[in] pdgSpin Spin.
     *  @param[in] electric_dipole_moment Electric dipole moment.
     *  @param[in] magnetic_charge Magnetic charge for monopoles.
     */
    inline G4ChargeState(G4double charge,
                         G4double magnetic_dipole_moment,
                         G4double pdgSpin, 
                         G4double electric_dipole_moment = 0.0,
                         G4double magnetic_charge = 0.0);

    /**
     * Copy constructor and assignment operator.
     */
    inline G4ChargeState( const G4ChargeState& right );
    inline G4ChargeState& operator = ( const G4ChargeState& right );

    /**
     * Default Destructor.
     */
    ~G4ChargeState() = default;

    /**
     * Revises the charge, pdgSpin, and optionally both moments and
     * magnetic charge.
     */
    void SetChargeSpinMoments(G4double charge,
                              G4double pdgSpin,  
                              G4double magnetic_dipole_moment= DBL_MAX,
                              G4double electric_dipole_moment= DBL_MAX,
                              G4double magnetic_charge= DBL_MAX );

    /**
     * Revises the charge (in units of the positron charge).
     */
    inline void SetCharge(G4double charge);
    inline G4double GetCharge() const;

    /**
     * Modifiers and accessors.
     */
    inline void SetPDGSpin(G4double spin);
    inline G4double GetPDGSpin() const;
    inline void SetSpin(G4double spin);
    inline G4double GetSpin() const;
    inline void SetMagneticDipoleMoment(G4double moment);
    inline G4double GetMagneticDipoleMoment() const;
    inline void SetElectricDipoleMoment(G4double moment);
    inline G4double ElectricDipoleMoment() const;
    inline void SetMagneticCharge(G4double charge);
    inline G4double MagneticCharge() const;

    /**
     * Auxiliary methods to set several properties at once.
     */
    inline void SetChargeMdm(G4double charge, G4double mag_dipole_moment);
    inline void SetChargeMdmSpin(G4double charge,
                                 G4double magnetic_dipole_moment,
                                 G4double pdgSpin);
    inline void SetChargeSpin(G4double charge,
                              G4double pdgSpin); 
    inline void SetChargeDipoleMoments(G4double charge,
                                       G4double magnetic_dipole_moment,
                                       G4double electric_dipole_moment);
    inline void SetChargesAndMoments(G4double charge,
                                     G4double magnetic_dipole_moment, 
                                     G4double electric_dipole_moment,
                                     G4double magnetic_charge );

   private:

     G4double fCharge;
     G4double fSpin;
     G4double fMagn_dipole;
     G4double fElec_dipole;
     G4double fMagneticCharge;  // for magnetic monopole
};

// Inline methods implementation

#include "G4ChargeState.icc"

#endif

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
// G4FieldTrack
//
// Class description:
//
// Data structure bringing together a magnetic track's state
// (position, momentum direction & modulus, energy, spin, ... ).
// Uses/abilities:
//  - does not maintain any relationship between its data (eg energy/momentum).
//  - for use in Runge-Kutta solver (in passing it the values right now).

// Author: John Apostolakis (CERN), 14.10.1996 - First version
// -------------------------------------------------------------------
#ifndef G4FIELDTRACK_HH
#define G4FIELDTRACK_HH

#include "G4ThreeVector.hh"
#include "G4ChargeState.hh"

/**
 * @brief G4FieldTrack defines a data structure bringing together a magnetic
 * track's state (position, momentum direction & modulus, energy, spin, etc. ).
 */

class G4FieldTrack
{
  public:

    /**
     * Constructor for G4FieldTrack.
     *  @param[in] pPosition Position in Cartesian coordinates.
     *  @param[in] LaboratoryTimeOfFlight Laboratory time of flight value.
     *  @param[in] pMomentumDirection Direction vector.
     *  @param[in] kineticEnergy Kinetic energy value.
     *  @param[in] restMass_c2 Mass at rest.
     *  @param[in] charge Charge.
     *  @param[in] polarization Polarisation vector.
     *  @param[in] magnetic_dipole_moment Magnetic dipole moment.
     *  @param[in] curve_length Length of curve.
     *  @param[in] PDGspin Spin.
     */
    G4FieldTrack( const G4ThreeVector& pPosition, 
                        G4double LaboratoryTimeOfFlight,
                  const G4ThreeVector& pMomentumDirection,
                        G4double kineticEnergy,
                        G4double restMass_c2,
                        G4double charge, 
                  const G4ThreeVector& polarization,
                        G4double magnetic_dipole_moment = 0.0,
                        G4double curve_length = 0.0,
                        G4double PDGspin = -1.0 );

    /**
     * Older constructor for G4FieldTrack, similar to above but missing charge.
     *  @param[in] pPosition Position in Cartesian coordinates.
     *  @param[in] pMomentumDirection Direction vector.
     *  @param[in] curve_length Length of curve.
     *  @param[in] kineticEnergy Kinetic energy value.
     *  @param[in] restMass_c2 Mass at rest.
     *  @param[in] velocity Velocity value - Not used.
     *  @param[in] LaboratoryTimeOfFlight Laboratory time of flight value.
     *  @param[in] ProperTimeOfFlight Proper time of flight value.
     *  @param[in] polarization Polarisation vector.
     *  @param[in] PDGspin Spin.
     */
    G4FieldTrack( const G4ThreeVector& pPosition, 
                  const G4ThreeVector& pMomentumDirection,
                        G4double       curve_length,
                        G4double       kineticEnergy,
                  const G4double       restMass_c2,
                        G4double       velocity,
                        G4double       LaboratoryTimeOfFlight = 0.0,
                        G4double       ProperTimeOfFlight = 0.0, 
                  const G4ThreeVector* pPolarization = nullptr,
                        G4double       PDGspin = -1.0 );

    /**
     * Empty init constructor.
     */
    G4FieldTrack( char );

    /**
     * Default Destructor.
     */
    ~G4FieldTrack() = default;

    /**
     * Copy constructor and assignment operator.
     */
    inline G4FieldTrack( const G4FieldTrack& pFieldTrack ); 
    inline G4FieldTrack& operator= ( const G4FieldTrack& rStVec );

    /**
     * Move constructor and move assignment operator.
     */
    inline G4FieldTrack(G4FieldTrack&& from) noexcept ;
    inline G4FieldTrack& operator=(G4FieldTrack&& from) noexcept ;

    /**
     * Streaming operator.
     */
    friend std::ostream& operator<<(std::ostream& os, const G4FieldTrack& SixVec);

    /**
     * Updates four-vectors for space/time and momentum/energy, also
     * resets the curve length.
     *  @param[in] pPosition Position in Cartesian coordinates.
     *  @param[in] LaboratoryTimeOfFlight Laboratory time of flight value.
     *  @param[in] pMomentumDirection Direction vector.
     *  @param[in] kineticEnergy Kinetic energy value.
     */
    inline void UpdateState( const G4ThreeVector& pPosition, 
                                   G4double LaboratoryTimeOfFlight,
                             const G4ThreeVector& pMomentumDirection,
                                   G4double kineticEnergy); 

    /**
     * Updates momentum, direction and kinetic energy.
     *  @param[in] kineticEnergy Kinetic energy value.
     *  @param[in] pMomentumDirection Direction vector.
     */
    inline void UpdateFourMomentum( G4double kineticEnergy, 
                                    const G4ThreeVector& momentumDirection ); 

    /**
     * Sets the charges and moments that are not given as DBL_MAX.
     *  @param[in] charge Charge value.
     *  @param[in] magnetic_dipole_moment M agnetic dipole moment.
     *  @param[in] electric_dipole_moment Electric dipole moment.
     *  @param[in] magnetic_charge Magnetic charge.
     */
    void SetChargeAndMoments(G4double charge, 
                             G4double magnetic_dipole_moment = DBL_MAX,
                             G4double electric_dipole_moment = DBL_MAX,
                             G4double magnetic_charge = DBL_MAX );

    /**
     * Setter and getter for PDG spin.
     */
    inline void SetPDGSpin(G4double pdgSpin);
    inline G4double GetPDGSpin();

    /**
     * Accessors.
     */
    inline G4ThreeVector GetMomentum() const;   
    inline G4ThreeVector GetPosition() const; 
    inline const G4ThreeVector& GetMomentumDir() const;
    inline G4ThreeVector GetMomentumDirection() const;
    inline G4double GetCurveLength() const;
    inline const G4ChargeState* GetChargeState() const;
    inline G4double GetLabTimeOfFlight() const;
    inline G4double GetProperTimeOfFlight() const;
    inline G4double GetKineticEnergy() const;
    inline G4double GetCharge() const;
    inline G4double GetRestMass() const;

    /**
     * Getter and setter for polarisation.
     */
    inline G4ThreeVector GetPolarization() const; 
    inline void SetPolarization( const G4ThreeVector& vecPol );

    /**
     * Setters for momentum. SetMomentumDir() does not change momentum
     * or Velocity Vector.
     */
    inline void SetMomentum(const G4ThreeVector& nMomDir);
    inline void SetMomentumDir(const G4ThreeVector& nMomDir);

    /**
     * Modifiers.
     */
    inline void SetPosition(const G4ThreeVector& nPos); 
    inline void SetRestMass(G4double Mass_c2);
    inline void SetCurveLength(G4double nCurve_s); // Distance along curve.
    inline void SetKineticEnergy(G4double nEnergy); // Does not modify momentum.
    inline void SetLabTimeOfFlight(G4double tofLab); 
    inline void SetProperTimeOfFlight(G4double tofProper);

    enum { ncompSVEC = 12 }; // Needed; should be used only for RK integration driver

    /**
     * Dumps/loads values to/from a provided array 'valArray'.
     */
    inline void DumpToArray(G4double valArr[ncompSVEC]) const; 
    void LoadFromArray(const G4double valArr[ncompSVEC],
                             G4int noVarsIntegrated);

    /**
     * More setters/getter foe spin, now obsolete.
     */
    inline void InitialiseSpin( const G4ThreeVector& vecPolarization );
    inline G4ThreeVector GetSpin() const;
    inline void SetSpin(const G4ThreeVector& vSpin);

  private:

    /**
     * Implementation method. Obsolete.
     */
    inline G4FieldTrack& SetCurvePnt(const G4ThreeVector& pPosition, 
                                     const G4ThreeVector& pMomentum,
                                           G4double       s_curve );
  private:

    G4double SixVector[6];
    G4double fDistanceAlongCurve;  // distance along curve of point
    G4double fKineticEnergy;
    G4double fRestMass_c2;
    G4double fLabTimeOfFlight;
    G4double fProperTimeOfFlight;
    G4ThreeVector fPolarization;
    G4ThreeVector fMomentumDir;
    G4ChargeState fChargeState;
}; 

#include "G4FieldTrack.icc"

#endif

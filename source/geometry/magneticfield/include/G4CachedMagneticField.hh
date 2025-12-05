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
// G4CachedMagneticField
//
// Class description:
//
// Caches Magnetic Field value, for field whose evaluation is expensive.

// Author: John Apostolakis (CERN), 20.07.2009.
// --------------------------------------------------------------------
#ifndef G4CACHED_MAGNETIC_FIELD_HH
#define G4CACHED_MAGNETIC_FIELD_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4MagneticField.hh"

/**
 * @brief G4CachedMagneticField is a specialisation of G4MagneticField and
 * is used to cache the Magnetic Field value, for fields whose evaluation is
 * expensive.
 */

class G4CachedMagneticField : public G4MagneticField
{
  public:

    /**
     * Constructor for G4CachedMagneticField.
     *  @param[in] pMagField Pointer to the original magnetic field.
     *  @param[in] distance Distance for field evaluation, within
     *             which the field does not change.
     */
    G4CachedMagneticField(G4MagneticField* pMagField, G4double distance);

    /**
     * Default Destructor.
     */
    ~G4CachedMagneticField() override = default;

    /**
     * Copy constructor and assignment operator.
     */
    G4CachedMagneticField(const G4CachedMagneticField& r);
    G4CachedMagneticField& operator = (const G4CachedMagneticField& p);

    /**
     * Returns the value of the field at the give 'Point'.
     *  @param[in] Point The given position time vector (x,y,z,t).
     *  @param[out] Bfield The returned field array.
     */
    void GetFieldValue( const G4double Point[4],
                              G4double* Bfield ) const override;
     
    /**
     * Getter and setter for the distance within which field is constant.
     */
    inline G4double GetConstDistance() const { return fDistanceConst; } 
    inline void SetConstDistance( G4double dist ) { fDistanceConst = dist;}

    /**
     * Accessors.
     */
    inline G4int GetCountCalls() const { return fCountCalls; }
    inline G4int GetCountEvaluations() const { return fCountEvaluations; } 

    /**
     * Resets counters.
     */
    inline void ClearCounts() { fCountCalls = 0; fCountEvaluations=0; }

    /**
     * Streams on standard output the values of counters.
     */
    void ReportStatistics();
    
    /**
     * Returns a pointer of an allocated clone of the field.
     */
    G4Field* Clone() const override;

  protected:

    mutable G4int fCountCalls = 0, fCountEvaluations = 0;  

  private:

    G4MagneticField* fpMagneticField = nullptr;

    /** When the field is evaluated within this distance it will not change. */
    G4double fDistanceConst;

    /** Caching state. */
    mutable G4ThreeVector fLastLocation;
    mutable G4ThreeVector fLastValue;
};

#endif

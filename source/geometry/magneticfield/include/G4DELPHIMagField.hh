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
// G4DELPHIMagField
//
// Class description:
//
// Class describing the DELPHI magnetic field. The axial symmetry
// field is mainly directed along Z axis. The function MagneticField(yTrack,B)
// calculates the magnetic induction vector B in point corresponding to
// yTrack according to parametrization given in:
//   P.Billoir, Precise tracking in a quasi-homogeneous magnetic field,
//              DELPHI 87-6 PROG 65, 1987.

// Author: Vladimir Grichine (CERN), 03.02.1997
// -------------------------------------------------------------------
#ifndef G4DELPHIMAGFIELD_HH
#define G4DELPHIMAGFIELD_HH

#include "G4MagneticField.hh"

/**
 * @brief describes the DELPHI magnetic field. The axial symmetry field is
 * mainly directed along Z axis. The function MagneticField(yTrack,B)
 * calculates the magnetic induction vector B in given point corresponding
 * according to parameterisation given in: P.Billoir, DELPHI 87-6 PROG 65, 1987.
 */

class G4DELPHIMagField : public G4MagneticField
{
  public:
                       
    /**
     * Default Constructor and Destructor.
     */
    G4DELPHIMagField() = default;
    ~G4DELPHIMagField() override = default;

    /**
     * Returns the field value on the given position 'yTrack'.
     *  @param[in] yTrack Time position array.
     *  @param[out] B The returned field array.
     */
    void GetFieldValue(const G4double yTrack[], G4double B[]) const override;

    /**
     * Returns a pointer to a new allocated clone of this object.
     */
    G4Field* Clone() const override;
};

#endif

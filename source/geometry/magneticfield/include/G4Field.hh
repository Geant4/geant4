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
// G4Field
//
// Class description:
//
// Abstract class for any kind of Field.
// It allows any kind of field (vector, scalar, tensor and any set of them)
// to be defined by implementing the inquiry function interface.
//
// The key method is  GetFieldValue( const G4double Point[4],
//                    *************        G4double* fieldArr ) 
// Given an input position/time vector 'Point', 
// this method must return the value of the field in "fieldArr".
//
// A field must also specify whether it changes a track's energy:
//                    DoesFieldChangeEnergy() 
//                    *********************
// A field must co-work with a corresponding Equation of Motion, to
// enable the integration of a particle's position, momentum and, optionally, 
// spin.  For this a field and its equation of motion must follow the
// same convention for the order of field components in the array "fieldArr"

// Author: John Apostolakis (CERN), 10.03.1997
// -------------------------------------------------------------------
#ifndef G4FIELD_HH
#define G4FIELD_HH

#include "G4Types.hh"
#include "G4FieldParameters.hh"
#include "globals.hh"

/**
 * @brief G4Field is the abstract class for any kind of field.
 * It allows any kind of field (vector, scalar, tensor and any set of them)
 * to be defined by implementing the inquiry function interface.
 * A field must co-work with a corresponding Equation of Motion, to
 * enable the integration of a particle's position, momentum and, optionally, 
 * spin. For this a field and its equation of motion must follow the same
 * convention for the order of field components.
 */

class G4Field
{
  public:

    /**
     * Constructor for G4Field.
     *  @param[in] gravityOn Flag to indicate if gravity is enabled or not.
     */
    G4Field(G4bool gravityOn = false);

    /**
     * Default virtual Destructor.
     */
    virtual ~G4Field() = default;

    /**
     * Copy constructor and assignment operator.
     */
    G4Field( const G4Field& p) = default;
    G4Field& operator = (const G4Field& p); 

    /**
     * Given the position time vector 'Point', returns the value of the
     * field in the array 'fieldArr'. Notes: 
     * 1) The 'Point' vector has the following structure:
     *      Point[0]  is  x  ( position, in Geant4 units )
     *      Point[1]  is  y
     *      Point[2]  is  z
     *      Point[3]  is  t  ( time, in Geant4 units )
     * 2) The convention for the components of the field array 'fieldArr'
     *    are determined by the type of field.
     *  @param[in] Point The position time vector.
     *  @param[out] fieldArr The field array in output.
     */
    virtual void GetFieldValue( const G4double Point[4],
                                      G4double* fieldArr ) const = 0;

    /**
     * Each type/class of field should respond the field does change energy.
     * For example:
     *  - an electric field     should return "true"
     *  - a pure magnetic field should return "false"
     */
    virtual G4bool DoesFieldChangeEnergy() const = 0;

    /**
     * Returns the field type-ID, "kUserFieldType".
     * This should be overriden in derived classes.
     */
    virtual G4FieldType GetFieldType() const { return kUserFieldType; }

    /**
     * Replies if the field includes gravity.
     *  @returns true if the field does include gravity.
     */
    inline G4bool IsGravityActive() const { return fGravityActive; }
        //  Does this field include gravity?

    /**
     * Sets the gravity flag.
     */
    inline void SetGravityActive(G4bool OnOffFlag) { fGravityActive = OnOffFlag; }
    
    /**
     * Interface method to implement cloning, needed by multi-threading.
     * Here issuing a fatal exception, as expecting this to be implemented
     * concretely in derived classes.
     */
    virtual G4Field* Clone() const;

  public:

    static constexpr G4int MAX_NUMBER_OF_COMPONENTS = 24;

  private:

    G4bool fGravityActive = false;
};
 
#endif

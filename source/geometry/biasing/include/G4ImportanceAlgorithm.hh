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
// G4ImportanceAlgorithm
//
// Class description:
//
// An implementation of a G4VImportanceAlgorithm (see description in
// G4VImportanceAlgorithm header).

// Author: Michael Dressel (CERN), 2002
// ----------------------------------------------------------------------
#ifndef G4IMPORTANCEALGORITHM_HH
#define G4IMPORTANCEALGORITHM_HH

#include "G4Threading.hh"
#include "G4VImportanceAlgorithm.hh"

/**
 * @brief G4ImportanceAlgorithm is a concrete implementation of a
 * G4VImportanceAlgorithm.
 */

class G4ImportanceAlgorithm : public G4VImportanceAlgorithm
{
  public:

    /**
     * Default Constructor and Destructor.
     */
    G4ImportanceAlgorithm() = default;
    ~G4ImportanceAlgorithm() override = default;

    /**
     * Calculates the number of tracks and their weight according to the
     * pre and post importance value and the weight of the mother track.
     *  @param[in] ipre "pre" importance value.
     *  @param[in] ipost "post" importance value.
     *  @param[in] init_w Initial weight value.
     *  @returns A struct containing the number of copies (including the
     *           mother track) to be produced and the weight of each track.
     */
    G4Nsplit_Weight Calculate(G4double ipre, 
                              G4double ipost, 
                              G4double init_w) const override;

  private:

    /**
     * Simple internal loggers.
     */
    void Error(const G4String& m) const;
    void Warning(const G4String& m) const;

  private:

    mutable G4bool fWarned = false;
};

#endif

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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLAllocationPool.hh"

#ifndef TransmissionChannel_hh
#define TransmissionChannel_hh 1
namespace G4INCL {
  class FinalState;

  class TransmissionChannel : public IChannel {
  public:
    TransmissionChannel(Nucleus *n, Particle *p);
    TransmissionChannel(Nucleus *n, Particle *p, const G4double TOut);
    TransmissionChannel(Nucleus *n, Particle *p, const G4double kOut, const G4double cosR);
    virtual ~TransmissionChannel();

    void fillFinalState(FinalState *fs);

  private:
    /** \brief Modify particle that leaves the nucleus.
     *
     * Modify the particle momentum and/or position when the particle leaves
     * the nucleus.
     */
    void particleLeaves();

    /** \brief Kinetic energy of the transmitted particle
     *
     * Calculate the kinetic energy of the particle outside the nucleus, if the
     * value has not been provided as a pre-calculated argument to the
     * constructor.
     */
    G4double initializeKineticEnergyOutside();

    Nucleus * const theNucleus;
    Particle * const theParticle;

    /// \brief True if refraction should be applied
    const G4bool refraction;

    /// \brief Momentum of the particle outside the nucleus
    const G4double pOutMag;

    /// \brief Kinetic energy of the particle outside the nucleus
    const G4double kineticEnergyOutside;

    /// \brief Cosine of the refraction angle
    const G4double cosRefractionAngle;

    INCL_DECLARE_ALLOCATION_POOL(TransmissionChannel)
  };
}
#endif // TransmissionChannel_hh

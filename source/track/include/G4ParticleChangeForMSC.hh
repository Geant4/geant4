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
// G4ParticleChangeForMSC
//
// Class description:
//
// Concrete "Particle Change" class for Multiple Scattering process.

// Author: Hisaya Kurashige, 23 March 1998  
// Revision: Vladimir Ivantchenko, 16 January 2004
// --------------------------------------------------------------------
#ifndef G4ParticleChangeForMSC_hh
#define G4ParticleChangeForMSC_hh 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include "G4ThreeVector.hh"
#include "G4VParticleChange.hh"

class G4DynamicParticle;

class G4ParticleChangeForMSC : public G4VParticleChange
{
  public:

    G4ParticleChangeForMSC();
      // Default constructor

    virtual ~G4ParticleChangeForMSC();
      // Destructor

  // --- the following methods are for updating G4Step -----
  // Return the pointer to the G4Step after updating the Step information
  // by using final state information of the track given by a physics
  // process

    virtual G4Step* UpdateStepForAlongStep(G4Step* Step);
    virtual G4Step* UpdateStepForPostStep(G4Step* Step);
      // A physics process gives the final state of the particle
      // based on information of G4Track (or equivalently the PreStepPoint)

    virtual void Initialize(const G4Track&);
      // Initialize all properties by using G4Track information

  // --- methods to keep information of the final state ---
  // IMPORTANT NOTE: despite the name, what the class stores through its
  //                 methods are the "FINAL" values of the Position,
  //                 Momentum, etc.

    void ProposeMomentumDirection(const G4ThreeVector& Pfinal);
    void ProposeMomentumDirection(G4double Px, G4double Py, G4double Pz);
    const G4ThreeVector* GetMomentumDirection() const;
    const G4ThreeVector* GetProposedMomentumDirection() const;
    void SetProposedMomentumDirection(const G4ThreeVector& Pfinal);
      // Get/Set theMomentumDirectionChange vector: it is the final momentum
      // direction

    const G4ThreeVector* GetPosition() const;
    void ProposePosition(const G4ThreeVector& finalPosition);
    const G4ThreeVector* GetProposedPosition() const;
    void SetProposedPosition(const G4ThreeVector& finalPosition);
      // Get/Set the final position of the current particle

  // --- Dump and debug methods ---

    virtual void DumpInfo() const;

    virtual G4bool CheckIt(const G4Track&);

  protected:

    G4ParticleChangeForMSC(const G4ParticleChangeForMSC& right);
    G4ParticleChangeForMSC& operator=(const G4ParticleChangeForMSC& right);
      // Hidden copy constructor and assignment operator

  private:

    G4ThreeVector theMomentumDirection;
      // It is the vector containing the final momentum direction
      // after the invoked process. The application of the change
      // of the momentum direction of the particle is not done here.
      // The responsibility to apply the change is up the entity
      // which invoked the process

    G4ThreeVector thePosition;
      // The changed (final) position of a given particle
};

#include "G4ParticleChangeForMSC.icc"

#endif

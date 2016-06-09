//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ParticleChangeForMSC.hh,v 1.6 2004/01/20 15:29:41 vnivanch Exp $
// GEANT4 tag $ $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// 
// Class Description
//   This class is special "Particle Change" for Multiple Scattering process
//
// ------------------------------------------------------------
//   Implemented for the new scheme                 23 Mar. 1998  H.Kurahige
//   Add Get/SetMomentumDirectionChange              6 Feb. 1999  H.Kurashige 
//   Update for model variant of msc                16 Jan  2004  V.Ivanchenko
//
// -------------------------------------------------------------
#ifndef G4ParticleChangeForMSC_h
#define G4ParticleChangeForMSC_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include "G4ThreeVector.hh"
class G4DynamicParticle;
#include "G4VParticleChange.hh"

class G4ParticleChangeForMSC: public G4VParticleChange
{ 
  public:
    // default constructor
    G4ParticleChangeForMSC();

    // destructor
    virtual ~G4ParticleChangeForMSC();

  protected:
    // hide copy constructor and assignment operaor as protected
    G4ParticleChangeForMSC(const G4ParticleChangeForMSC &right);
    G4ParticleChangeForMSC & operator=(const G4ParticleChangeForMSC &right);


  public: // with description
    // ----------------------------------------------------
    // --- the following methods are for updating G4Step -----
    // Return the pointer to the G4Step after updating the Step information
    // by using final state information of the track given by a physics
    // process
    virtual G4Step* UpdateStepForAlongStep(G4Step* Step);
    virtual G4Step* UpdateStepForPostStep(G4Step* Step);
    // A physics process gives the final state of the particle
    // based on information of G4Track (or equivalently the PreStepPoint)

    virtual void Initialize(const G4Track&);
    // Initialize all propoerties by using G4Track information

    // ----------------------------------------------------
    //--- methods to keep information of the final state--
    //  IMPORTANT NOTE: Although the name of the class and methods are
    //   "Change", what it stores (and returns in get) are the "FINAL"
    //   values of the Position, Momentum, etc.

    void SetMomentumChange(const G4ThreeVector& Pfinal);
    void SetMomentumChange(G4double Px, G4double Py, G4double Pz);
    const G4ThreeVector* GetProposedMomentumDirection() const;
    void SetProposedMomentumDirection(const G4ThreeVector& Pfinal);
    // Get/Set theMomentumDirectionChange vector: it is the final momentum direction.

    const G4ThreeVector* GetPositionChange() const;
    void SetPositionChange(const G4ThreeVector& finalPosition);
    const G4ThreeVector* GetProposedPosition() const;
    void SetProposedPosition(const G4ThreeVector& finalPosition);
    //  Get/Set the final position of the current particle.

  public:
    virtual void DumpInfo() const;
    // for Debug
    virtual G4bool CheckIt(const G4Track&);

  private:
    G4ThreeVector theMomentumDirection;
    //  It is the vector containing the final momentum direction
    //  after the invoked process. The application of the change
    //  of the momentum direction of the particle is not Done here.
    //  The responsibility to apply the change is up the entity
    //  which invoked the process.

    G4ThreeVector thePosition;
    //  The changed (final) position of a given particle.

};

inline
 void G4ParticleChangeForMSC::SetMomentumChange(const G4ThreeVector& P)
{
  theMomentumDirection = P;
}

inline
 void G4ParticleChangeForMSC::SetProposedMomentumDirection(const G4ThreeVector& P)
{
  theMomentumDirection = P;
}

inline
 void G4ParticleChangeForMSC::SetMomentumChange(G4double Px, G4double Py, G4double Pz)
{
  theMomentumDirection.setX(Px);
  theMomentumDirection.setY(Py);
  theMomentumDirection.setZ(Pz);
}

inline
 const G4ThreeVector* G4ParticleChangeForMSC::GetProposedMomentumDirection() const
{
  return &theMomentumDirection;
}

inline
 void G4ParticleChangeForMSC::SetProposedPosition(const G4ThreeVector& P)
{
  thePosition = P;
}

inline
 const G4ThreeVector* G4ParticleChangeForMSC::GetProposedPosition() const
{
  return &thePosition;
}

inline
 void G4ParticleChangeForMSC::SetPositionChange(const G4ThreeVector& P)
{
  thePosition = P;
}

inline
 const G4ThreeVector* G4ParticleChangeForMSC::GetPositionChange() const
{
  return &thePosition;
}

inline void G4ParticleChangeForMSC::Initialize(const G4Track& track)
{
  theStatusChange = track.GetTrackStatus();
  theMomentumDirection = track.GetMomentumDirection();
  thePosition = track.GetPosition();
}

#endif

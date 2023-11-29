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
//

#ifndef G4VSplitableHadron_h
#define G4VSplitableHadron_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4VSplitableHadron----------------
//             by Gunter Folger, June 1998.
//       class storing an interacting particle. Used by Parton String Models.
// ------------------------------------------------------------

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ReactionProduct.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

class G4Nucleon;
class G4Parton;
class G4VKineticNucleon;

#include <vector>

class G4VSplitableHadron 
{
  public:
    G4VSplitableHadron();
    G4VSplitableHadron(const G4ReactionProduct & aPrimary);
    G4VSplitableHadron(const G4Nucleon & aNucleon);
    G4VSplitableHadron(const G4VKineticNucleon * aNucleon);

    virtual ~G4VSplitableHadron();

    G4bool operator==(const G4VSplitableHadron &right) const;
    G4bool operator!=(const G4VSplitableHadron &right) const;

    void Set4Momentum(const G4LorentzVector &a4Momentum);
    const G4LorentzVector & Get4Momentum() const;

    void SetDefinition(const G4ParticleDefinition *aDefinition);
    const G4ParticleDefinition * GetDefinition() const;
      
    void IncrementCollisionCount(G4int aCount);
    void SetCollisionCount(G4int aCount);

    void SetTimeOfCreation(G4double aTime);
    G4double GetTimeOfCreation();

    void SetPosition(const G4ThreeVector &aPosition);
    const G4ThreeVector & GetPosition() const;

    void SetStatus(const G4int aStatus);
    G4int GetStatus();

    virtual void SplitUp() = 0;
    virtual void SetFirstParton(G4int PDGcode) = 0;
    virtual void SetSecondParton(G4int PDGcode)= 0;
    virtual G4Parton * GetNextParton() = 0 ;
    virtual G4Parton * GetNextAntiParton() = 0 ;
    G4bool IsSplit() { return isSplit;}

    G4int GetSoftCollisionCount();

    void Splitting() {isSplit = true;}

  private:
    G4VSplitableHadron(const G4VSplitableHadron &right);
    const G4VSplitableHadron & operator=(const G4VSplitableHadron &right);

    const G4ParticleDefinition *theDefinition;

    G4LorentzVector the4Momentum;

    G4double TimeOfCreation;
    G4ThreeVector thePosition;
    G4int theCollisionCount;

    G4int  curStatus;
    G4bool isSplit;
};

inline G4int G4VSplitableHadron::GetSoftCollisionCount()
{
  return theCollisionCount;
}

inline void G4VSplitableHadron::SetCollisionCount(G4int aCount)
{
  theCollisionCount = aCount;
}

inline void G4VSplitableHadron::Set4Momentum(const G4LorentzVector &a4Momentum)
{
  the4Momentum=a4Momentum;
}

inline const G4LorentzVector & G4VSplitableHadron::Get4Momentum() const
{
  return the4Momentum;
}

inline void G4VSplitableHadron::SetDefinition(const G4ParticleDefinition *aDefinition)
{
  theDefinition=aDefinition;
}

inline const G4ParticleDefinition * G4VSplitableHadron::GetDefinition() const
{
  return theDefinition;
}

inline void G4VSplitableHadron::IncrementCollisionCount(G4int aCount)
{
  theCollisionCount += aCount;
}

inline void G4VSplitableHadron::SetTimeOfCreation(G4double aTime)
{
  TimeOfCreation=aTime;
}

inline G4double G4VSplitableHadron::GetTimeOfCreation()
{
  return TimeOfCreation; 
}

inline void G4VSplitableHadron::SetPosition(const G4ThreeVector &aPosition)
{
  thePosition=aPosition;
}

inline const G4ThreeVector & G4VSplitableHadron::GetPosition() const
{
  return thePosition;
}

inline void G4VSplitableHadron::SetStatus(G4int aStatus)
{
  curStatus=aStatus;
}

inline G4int G4VSplitableHadron::GetStatus()
{
  return curStatus; 
}

#endif


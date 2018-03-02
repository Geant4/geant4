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
// $Id: G4Cerenkov.hh 108508 2018-02-15 15:54:35Z gcosmo $
//
// 
////////////////////////////////////////////////////////////////////////
// Cerenkov Radiation Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4Cerenkov.hh
// Description:	Discrete Process - Generation of Cerenkov Photons
// Version:     2.0
// Created:     1996-02-21
// Author:      Juliet Armstrong
// Updated:     2007-09-30 change inheritance to G4VDiscreteProcess
//              2005-07-28 add G4ProcessType to constructor
//              1999-10-29 add method and class descriptors
//              1997-04-09 by Peter Gumplinger
//              > G4MaterialPropertiesTable; new physics/tracking scheme
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef G4Cerenkov_h
#define G4Cerenkov_h 1

/////////////
// Includes
/////////////

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"

// Class Description:
// Discrete Process -- Generation of Cerenkov Photons.
// Class inherits publicly from G4VDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class G4Cerenkov : public G4VProcess
{

public:

  ////////////////////////////////
  // Constructors and Destructor
  ////////////////////////////////

  explicit G4Cerenkov(const G4String& processName = "Cerenkov", 
             G4ProcessType type = fElectromagnetic);
  ~G4Cerenkov();

  explicit G4Cerenkov(const G4Cerenkov &right);

private:

  //////////////
  // Operators
  //////////////

  G4Cerenkov& operator=(const G4Cerenkov &right) = delete;

public:

  ////////////
  // Methods
  ////////////

  G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;
  // Returns true -> 'is applicable', for all charged particles
  // except short-lived particles.

  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType) override;
  // Build table at a right time

  G4double GetMeanFreePath(const G4Track& aTrack, 
                           G4double, G4ForceCondition* );
  // Returns the discrete step limit and sets the 'StronglyForced'
  // condition for the DoIt to be invoked at every step.

  G4double PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                                G4double ,
                                                G4ForceCondition* ) override;
  // Returns the discrete step limit and sets the 'StronglyForced'
  // condition for the DoIt to be invoked at every step.

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
                                  const G4Step&  aStep) override;
  // This is the method implementing the Cerenkov process.

  //  no operation in  AtRestDoIt and  AlongStepDoIt
  virtual G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
                                                         G4double  ,
                                                         G4double  ,
                                                         G4double& ,
                                                         G4GPILSelection*
                                                ) override { return -1.0; };

  virtual G4double AtRestGetPhysicalInteractionLength(const G4Track& ,
                                                      G4ForceCondition*
                                                 ) override { return -1.0; };

  //  no operation in  AtRestDoIt and  AlongStepDoIt
  virtual G4VParticleChange* AtRestDoIt(const G4Track& , const G4Step& ) 
                             override {return nullptr;};

  virtual G4VParticleChange* AlongStepDoIt(const G4Track& , const G4Step&) 
                             override {return nullptr;};

  void SetTrackSecondariesFirst(const G4bool state);
  // If set, the primary particle tracking is interrupted and any 
  // produced Cerenkov photons are tracked next. When all have 
  // been tracked, the tracking of the primary resumes.

  G4bool GetTrackSecondariesFirst() const;
  // Returns the boolean flag for tracking secondaries first.

  void SetMaxBetaChangePerStep(const G4double d);
  // Set the maximum allowed change in beta = v/c in % (perCent)
  // per step.

  G4double GetMaxBetaChangePerStep() const;
  // Returns the maximum allowed change in beta = v/c in % (perCent)

  void SetMaxNumPhotonsPerStep(const G4int NumPhotons);
  // Set the maximum number of Cerenkov photons allowed to be 
  // generated during a tracking step. This is an average ONLY; 
  // the actual number will vary around this average. If invoked, 
  // the maximum photon stack will roughly be of the size set.
  // If not called, the step is not limited by the number of 
  // photons generated.

  G4int GetMaxNumPhotonsPerStep() const;
  // Returns the maximum number of Cerenkov photons allowed to be
  // generated during a tracking step.

  void SetStackPhotons(const G4bool );
  // Call by the user to set the flag for stacking the scint. photons

  G4bool GetStackPhotons() const;
  // Return the boolean for whether or not the scint. photons are stacked

  G4int GetNumPhotons() const;
  // Returns the current number of scint. photons (after PostStepDoIt)

  G4PhysicsTable* GetPhysicsTable() const;
  // Returns the address of the physics table.

  void DumpPhysicsTable() const;
  // Prints the physics table.

private:

  void BuildThePhysicsTable();

  /////////////////////
  // Helper Functions
  /////////////////////

  G4double GetAverageNumberOfPhotons(const G4double charge,
                                     const G4double beta,
                                     const G4Material *aMaterial,
                                     G4MaterialPropertyVector* Rindex) const;

  ///////////////////////
  // Class Data Members
  ///////////////////////

protected:

  G4PhysicsTable* thePhysicsTable;
  //  A Physics Table can be either a cross-sections table or
  //  an energy table (or can be used for other specific
  //  purposes).

private:

  G4bool fTrackSecondariesFirst;
  G4double fMaxBetaChange;
  G4int  fMaxPhotons;

  G4bool fStackingFlag;

  G4int fNumPhotons;
};

  ////////////////////
  // Inline methods
  ////////////////////

inline
G4bool G4Cerenkov::GetTrackSecondariesFirst() const
{
  return fTrackSecondariesFirst;
}

inline
G4double G4Cerenkov::GetMaxBetaChangePerStep() const
{
  return fMaxBetaChange;
}

inline
G4int G4Cerenkov::GetMaxNumPhotonsPerStep() const
{
  return fMaxPhotons;
}

inline
void G4Cerenkov::SetStackPhotons(const G4bool stackingFlag)
{
        fStackingFlag = stackingFlag;
}

inline
G4bool G4Cerenkov::GetStackPhotons() const
{
        return fStackingFlag;
}

inline
G4int G4Cerenkov::GetNumPhotons() const
{
        return fNumPhotons;
}

inline
G4PhysicsTable* G4Cerenkov::GetPhysicsTable() const
{
  return thePhysicsTable;
}

#endif /* G4Cerenkov_h */

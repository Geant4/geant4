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
// $Id: G4VAtomDeexcitation.hh 98985 2016-08-29 07:01:17Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4VAtomDeexcitation
//
// Author:        Alfonso Mantero & Vladimir Ivanchenko
//
// Creation date: 30.06.2009
//
// Modifications:
// 15 Mar 2011  ALF   stripped G4AtomicShellEnumerator to its own file
//
// Class Description:
//
// Abstract interface to energy loss models

// -------------------------------------------------------------------
//

#ifndef G4VAtomDeexcitation_h
#define G4VAtomDeexcitation_h 1

#include "globals.hh"
#include "G4EmParameters.hh"
#include "G4AtomicShell.hh"
#include "G4AtomicShellEnumerator.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Track.hh"
#include <vector>

class G4ParticleDefinition;
class G4DynamicParticle;
class G4MaterialCutsCouple;

class G4VAtomDeexcitation {
public:

  explicit G4VAtomDeexcitation(const G4String& modname = "Deexcitation");

  virtual ~G4VAtomDeexcitation();

  //========== initialization ==========

  // Overall initialisation before new run
  void InitialiseAtomicDeexcitation();

  // Initialisation of deexcitation at the beginning of run 
  virtual void InitialiseForNewRun() = 0;

  // Initialisation for a concrete atom 
  // May be called at run time 
  virtual void InitialiseForExtraAtom(G4int Z) = 0;

  void SetDeexcitationActiveRegion(const G4String& rname, 
                                   G4bool valDeexcitation,
                                   G4bool valAuger,
                                   G4bool valPIXE);

  // Activation of deexcitation
  inline void SetFluo(G4bool);
  inline G4bool IsFluoActive() const;

  // Activation of Auger electron production
  inline void SetAuger(G4bool);
  inline G4bool IsAugerActive() const;

  // Activation of Auger cascade
  inline void SetAugerCascade(G4bool);
  inline G4bool IsAugerCascadeActive() const;

  // Activation of PIXE simulation
  inline void SetPIXE(G4bool);
  inline G4bool IsPIXEActive() const;

  // Deexcitation model name
  inline const G4String& GetName() const;

  // Access to the list of atoms active for deexcitation
  inline const std::vector<G4bool>& GetListOfActiveAtoms() const;

  // Verbosity level
  inline void SetVerboseLevel(G4int);
  inline G4int GetVerboseLevel() const;

  //========== Run time methods ==========

  // Check if deexcitation is active for a given geometry volume
  inline G4bool CheckDeexcitationActiveRegion(G4int coupleIndex);
  inline G4bool CheckAugerActiveRegion(G4int coupleIndex);

  // Get atomic shell by shell index, used by discrete processes 
  // (for example, photoelectric), when shell vacancy sampled by the model
  virtual 
  const G4AtomicShell* GetAtomicShell(G4int Z, 
                                      G4AtomicShellEnumerator shell) = 0;

  // generation of deexcitation for given atom and shell vacancy
  inline void GenerateParticles(std::vector<G4DynamicParticle*>* secVect,  
                                const G4AtomicShell*, 
                                G4int Z,
                                G4int coupleIndex);

  // generation of deexcitation for given atom and shell vacancy
  virtual void GenerateParticles(std::vector<G4DynamicParticle*>* secVect,  
                                 const G4AtomicShell*, 
                                 G4int Z,
                                 G4double gammaCut,
                                 G4double eCut) = 0;

  // access or compute PIXE cross section 
  virtual G4double 
  GetShellIonisationCrossSectionPerAtom(const G4ParticleDefinition*, 
                                        G4int Z, 
                                        G4AtomicShellEnumerator shell,
                                        G4double kinE,
                                        const G4Material* mat = nullptr) = 0;

  // access or compute PIXE cross section 
  virtual G4double 
  ComputeShellIonisationCrossSectionPerAtom(
                                        const G4ParticleDefinition*, 
                                        G4int Z, 
                                        G4AtomicShellEnumerator shell,
                                        G4double kinE,
                                        const G4Material* mat = nullptr) = 0;

  // Sampling of PIXE for ionisation processes
  void AlongStepDeexcitation(std::vector<G4Track*>& tracks,  
                             const G4Step& step, 
                             G4double& eLoss,
                             G4int coupleIndex);

private:

  // copy constructor and hide assignment operator
  G4VAtomDeexcitation(G4VAtomDeexcitation &) = delete;
  G4VAtomDeexcitation & operator=(const G4VAtomDeexcitation &right) = delete;

  G4EmParameters*             theParameters;  
  const G4ParticleDefinition* gamma;

  G4ProductionCutsTable* theCoupleTable;
  G4int    verbose;
  G4String name;

  G4bool   isActive;
  G4bool   flagAuger;
  G4bool   flagAugerCascade;
  G4bool   flagPIXE;
  G4bool   ignoreCuts;

  G4bool   isActiveLocked;
  G4bool   isAugerLocked;
  G4bool   isAugerCascadeLocked;
  G4bool   isPIXELocked;

  std::vector<G4bool>   activeZ;
  std::vector<G4bool>   activeDeexcitationMedia;
  std::vector<G4bool>   activeAugerMedia;
  std::vector<G4bool>   activePIXEMedia;
  std::vector<G4String> activeRegions;
  std::vector<G4bool>   deRegions;
  std::vector<G4bool>   AugerRegions;
  std::vector<G4bool>   PIXERegions;
  std::vector<G4DynamicParticle*> vdyn;

  static G4int pixeIDg;
  static G4int pixeIDe;
};

inline void G4VAtomDeexcitation::SetFluo(G4bool val)
{
  if(!isActiveLocked) { isActive = val; isActiveLocked = true; }
}

inline G4bool G4VAtomDeexcitation::IsFluoActive() const
{
  return isActive;
}

inline void G4VAtomDeexcitation::SetAuger(G4bool val)
{
  if(!isAugerLocked) { flagAuger = val; isAugerLocked = true; }
}

inline G4bool G4VAtomDeexcitation::IsAugerActive() const
{
  return flagAuger;
}

inline void G4VAtomDeexcitation::SetAugerCascade(G4bool val)
{
  if(!isAugerCascadeLocked) { flagAugerCascade = val;  isAugerCascadeLocked = true; }
}

inline G4bool G4VAtomDeexcitation::IsAugerCascadeActive() const
{
  return flagAugerCascade;
}

inline void G4VAtomDeexcitation::SetPIXE(G4bool val)
{
  if(!isPIXELocked) { flagPIXE = val;  isPIXELocked = true; }
}

inline G4bool G4VAtomDeexcitation::IsPIXEActive() const
{
  return flagPIXE;
}

inline const G4String& G4VAtomDeexcitation::GetName() const
{
  return name;
}

inline const std::vector<G4bool>& 
G4VAtomDeexcitation::GetListOfActiveAtoms() const
{
  return activeZ;
}

inline void G4VAtomDeexcitation::SetVerboseLevel(G4int val)
{
  verbose = val;
}

inline G4int G4VAtomDeexcitation::GetVerboseLevel() const
{
  return verbose;
}

inline G4bool 
G4VAtomDeexcitation::CheckDeexcitationActiveRegion(G4int coupleIndex)
{
  return (activeDeexcitationMedia[coupleIndex]);
}

inline G4bool 
G4VAtomDeexcitation::CheckAugerActiveRegion(G4int coupleIndex)
{
  return (activeAugerMedia[coupleIndex]);
}

inline void 
G4VAtomDeexcitation::GenerateParticles(std::vector<G4DynamicParticle*>* v,  
                                       const G4AtomicShell* as, 
                                       G4int Z,
                                       G4int idx)
{
  G4double gCut = DBL_MAX;
  if(ignoreCuts) {
    gCut = 0.0;
  } else if (theCoupleTable) {
    gCut = (*(theCoupleTable->GetEnergyCutsVector(0)))[idx];
  }
  if(gCut < as->BindingEnergy()) {
    G4double eCut = DBL_MAX;
    if(CheckAugerActiveRegion(idx)) {
      if(ignoreCuts) {
        eCut = 0.0;
      } else if (theCoupleTable) {
        eCut = (*(theCoupleTable->GetEnergyCutsVector(1)))[idx];
      }
    }
    GenerateParticles(v, as, Z, gCut, eCut);
  }
}

#endif


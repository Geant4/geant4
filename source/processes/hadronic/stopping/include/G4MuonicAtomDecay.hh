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
#ifndef G4MuonicAtomDecay_h
#define G4MuonicAtomDecay_h 1

//
//
// $Id: G4MuonicAtomDecay.hh 96314 2016-04-06 07:21:51Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: 
//
//  20170522 K L Genser first implementation somewhat based on
//  G4RadioactiveDecay and code by V.Ivantchenko & K.Lynch
//
// ------------------------------------------------------------

#include "G4VRestDiscreteProcess.hh"

class G4MuonicAtom;
class G4Ions;
class G4HadronicInteraction;
class G4HadFinalState;

class G4MuonicAtomDecay : public G4VRestDiscreteProcess
{

  // Class Description:
  //
  // MuonicAtom process in which Muon is captured by the nucleus or decays in orbit

public:  // with description

  explicit G4MuonicAtomDecay(G4HadronicInteraction* hiptr=nullptr,
                             // G4HadronicInteraction* diptr=nullptr,
                             const G4String& processName="MuonicAtomDecay");

  virtual ~G4MuonicAtomDecay();

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  //  virtual void PreparePhysicsTable(const G4ParticleDefinition&);

  //  virtual void BuildPhysicsTable(const G4ParticleDefinition&);

  inline void  SetVerboseLevel(G4int value) {verboseLevel = value;}
  // Sets the VerboseLevel

  inline G4int GetVerboseLevel() const {return verboseLevel;}
  // Returns the VerboseLevel

  virtual void ProcessDescription(std::ostream& outFile) const;

  virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                        G4double   previousStepSize,
                                                        G4ForceCondition* condition);

  virtual G4double AtRestGetPhysicalInteractionLength(const G4Track&,
                                                      G4ForceCondition* );
      
  inline virtual G4VParticleChange* AtRestDoIt(const G4Track& theTrack,
                                       const G4Step& theStep)
  {return DecayIt(theTrack, theStep);}
  
  inline virtual G4VParticleChange* PostStepDoIt(const G4Track& theTrack,
                                         const G4Step& theStep)
  {return DecayIt(theTrack, theStep);}

  virtual G4double GetMeanFreePath(const G4Track& aTrack,
                                   G4double   previousStepSize,
                                   G4ForceCondition* condition
                                   );
  //  Calculates from the macroscopic cross section a mean
  //  free path, the value is returned in units of distance.

  virtual G4double GetMeanLifeTime(const G4Track& aTrack,G4ForceCondition* condition);
  //  Calculates the mean life-time (i.e. for decays) of the
  //  particle at rest due to the occurence of the given process,
  //  or converts the probability of interaction (i.e. for
  //  annihilation) into the life-time of the particle for the
  //  occurence of the given process.


private:
  // hide copy constructor and assignment operator as private 
  G4MuonicAtomDecay(const G4MuonicAtomDecay &right);
  G4MuonicAtomDecay& operator=(const G4MuonicAtomDecay &right);

  G4VParticleChange* DecayIt(const G4Track& theTrack,
                             const G4Step&  theStep);  

  void FillResult(G4HadFinalState * aR, const G4Track & aT);

  void DumpState(const G4Track&, const G4String&, G4ExceptionDescription&);

  G4ParticleChange theTotalResult;

  const G4double fMuMass;
  
  G4HadronicInteraction* cmptr; // the capture model owned by its registry
  //  G4HadronicInteraction* dmptr; // the DIO model owned by its registry

  G4int verboseLevel;


};
#endif

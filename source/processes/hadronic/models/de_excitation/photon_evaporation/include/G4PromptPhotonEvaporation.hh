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
// $Id$
//
// -------------------------------------------------------------------
//
//      GEANT4 class file
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4PromptPhotonEvaporation
//
//      Author:        Vladimir Ivantchenko
//
//      Creation date: 20 December 2011
//
//Modifications:
//
// 
// -------------------------------------------------------------------

#ifndef G4PROMPTPHOTONEVAPORATION_HH
#define G4PROMPTPHOTONEVAPORATION_HH 1

#include "globals.hh"
#include "G4VEvaporationChannel.hh"

class G4Fragment;
class G4NuclearLevelManager;
class G4NuclearLevelStore;

class G4PromptPhotonEvaporation : public G4VEvaporationChannel {

public:

  G4PromptPhotonEvaporation();

  virtual ~G4PromptPhotonEvaporation();

  virtual G4double GetEmissionProbability(G4Fragment* theNucleus);

  // one photon emission
  virtual G4Fragment* EmittedFragment(G4Fragment* theNucleus);

  // full fragment de-excitation by photon emission
  virtual G4FragmentVector* BreakUpFragment(G4Fragment* theNucleus);

  // obsolete method
  virtual G4FragmentVector* BreakUp(const G4Fragment& theNucleus);

  inline void SetVerboseLevel(G4int verbose);

  inline void SetICM(G4bool);

  inline void RDMForced (G4bool);
  
  inline void SetMaxHalfLife(G4double);
 
private:

  G4PromptPhotonEvaporation(const G4PromptPhotonEvaporation & right);
  const G4PromptPhotonEvaporation & operator = (const G4PromptPhotonEvaporation & right);

  G4int    fVerbose;
  G4bool   fICM;
  G4bool   fRDM;
  G4double fMaxHalfTime;
  G4double fEmissionProbability;

  G4NuclearLevelManager* levelManager;
  G4NuclearLevelStore* fNuclearLevelStore;

  const G4Fragment* nucleus;

  G4int    theZ;
  G4int    theA;
  G4double fEnergyFermi;
  G4double fExcEnergyMax;

  G4double gammaE;

};

inline void G4PromptPhotonEvaporation::SetVerboseLevel(G4int verbose)
{
  fVerbose = verbose;
}

inline void G4PromptPhotonEvaporation::SetICM(G4bool val)
{
  fICM = val;
}

inline void G4PromptPhotonEvaporation::RDMForced (G4bool val)
{
  fRDM = val;
}
  
inline void G4PromptPhotonEvaporation::SetMaxHalfLife(G4double val)
{
  fMaxHalfTime = val;
}

#endif

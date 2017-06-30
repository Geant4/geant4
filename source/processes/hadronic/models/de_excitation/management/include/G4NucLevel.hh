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
// $Id: G4NucLevel.hh 88516 2015-02-25 11:00:16Z vnivanch $
//
// -------------------------------------------------------------------
//
//      GEANT4 header file 
//
//      File name:     G4NucLevel
//
//      Author:        V.Ivanchenko
// 
//      Creation date: 4 January 2012
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//  Container class keeping information about gamma transition
//  and for a given nuclear level
//

#ifndef G4NUCLEVEL_HH
#define G4NUCLEVEL_HH 1

#include "globals.hh"
#include <vector>

class G4NucLevel 
{
public:

  explicit G4NucLevel(size_t ntrans, G4double  tgamma,
	              const std::vector<G4int>&   vTrans,
	              const std::vector<G4float>& wLevelGamma,
	              const std::vector<G4float>& wGamma,
	              const std::vector<G4float>& vRatio,
                      const std::vector<const std::vector<G4float>*>& wShell);

  ~G4NucLevel();

  inline size_t NumberOfTransitions() const;

  inline size_t FinalExcitationIndex(size_t idx) const;

  inline G4int TransitionType(size_t idx) const;

  inline G4double GetTimeGamma() const;

  inline G4float MixingRatio(size_t idx) const;

  inline G4float GammaProbability(size_t idx) const;

  inline G4float GammaCumProbability(size_t idx) const;

  inline G4float MultipolarityRatio(size_t idx) const;

  inline size_t SampleGammaTransition(G4double rndm) const;

  inline G4int SampleShell(size_t idx, G4double rndm) const;

  inline const std::vector<G4float>* ShellProbabilty(size_t idx) const;

private:  

#ifdef G4VERBOSE
  void PrintError(size_t idx, const G4String&) const;  
#endif
  G4NucLevel(const G4NucLevel &right) = delete;
  G4bool operator==(const G4NucLevel &right) const = delete;
  G4bool operator!=(const G4NucLevel &right) const = delete;
  G4bool operator<(const G4NucLevel &right) const = delete;
  const G4NucLevel& operator=(const G4NucLevel &right) = delete;

  size_t   length;
  G4double fTimeGamma;
  
  std::vector<G4int>    fTrans;
  std::vector<G4float>  fGammaCumProbability;
  std::vector<G4float>  fGammaProbability;
  std::vector<G4float>  fMpRatio;
  std::vector<const std::vector<G4float>*> fShellProbability;
};

inline size_t G4NucLevel::NumberOfTransitions() const
{
  return length;
}

inline size_t G4NucLevel::FinalExcitationIndex(size_t idx) const
{
#ifdef G4VERBOSE
  if(idx >= length) { PrintError(idx, "FinalExcitationEnergy"); }
#endif
  return (size_t)(fTrans[idx]/10000);
}

inline G4int G4NucLevel::TransitionType(size_t idx) const
{
#ifdef G4VERBOSE
  if(idx >= length) { PrintError(idx, "TransitionType"); }
#endif
  return fTrans[idx]%10000;
}

inline G4double G4NucLevel::GetTimeGamma() const
{
  return fTimeGamma;
}

inline G4float G4NucLevel::MixingRatio(size_t idx) const
{
#ifdef G4VERBOSE
  if(idx >= length) { PrintError(idx, "MixingRatio"); }
#endif
  return fMpRatio[idx];
}

inline G4float G4NucLevel::GammaProbability(size_t idx) const
{
#ifdef G4VERBOSE
  if(idx >= length) { PrintError(idx, "GammaProbability"); }
#endif
  return fGammaProbability[idx];
}

inline G4float G4NucLevel::GammaCumProbability(size_t idx) const
{
#ifdef G4VERBOSE
  if(idx >= length) { PrintError(idx, "GammaCumProbability"); }
#endif
  return fGammaCumProbability[idx];
}

inline G4float G4NucLevel::MultipolarityRatio(size_t idx) const
{
#ifdef G4VERBOSE
  if(idx >= length) { PrintError(idx, "GammaProbability"); }
#endif
  return fMpRatio[idx];
}

inline size_t G4NucLevel::SampleGammaTransition(G4double rndm) const
{
  G4float x = (G4float)rndm;
  size_t idx = 0;
  for(; idx<length; ++idx) { 
    if(x <= fGammaCumProbability[idx]) { break; } 
  }
  return idx;
}

inline G4int G4NucLevel::SampleShell(size_t idx, G4double rndm) const
{
#ifdef G4VERBOSE
  if(idx >= length) { PrintError(idx, "SampleShell"); }
#endif
  const std::vector<G4float>* prob = fShellProbability[idx];
  G4int i(-1);
  if(prob) {
    G4int nn = prob->size();
    G4float x = (G4float)rndm;
    for(i=0; i<nn; ++i) { if(x <= (*prob)[i]) { break; } }
  } 
  return i;
}

inline const std::vector<G4float>* 
G4NucLevel::ShellProbabilty(size_t idx) const
{
#ifdef G4VERBOSE
  if(idx >= length) { PrintError(idx, "ShellProbability"); }
#endif
  return fShellProbability[idx];
}

#endif






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
#include <iostream>

class G4NucLevel 
{
public:

  explicit G4NucLevel(std::size_t ntrans, G4double  tgamma,
	              const std::vector<G4int>&   vTrans,
	              const std::vector<G4float>& wLevelGamma,
	              const std::vector<G4float>& wGamma,
	              const std::vector<G4float>& vRatio,
                      const std::vector<const std::vector<G4float>*>& wShell);

  ~G4NucLevel();

  inline std::size_t NumberOfTransitions() const;

  inline std::size_t FinalExcitationIndex(std::size_t idx) const;

  inline G4int TransitionType(std::size_t idx) const;

  inline G4double GetTimeGamma() const;

  inline G4float GammaProbability(std::size_t idx) const;

  inline G4float GammaCumProbability(std::size_t idx) const;

  inline G4float MultipolarityRatio(std::size_t idx) const;

  inline std::size_t SampleGammaTransition(G4double rndm) const;

  inline G4int SampleShell(std::size_t idx, G4double rndm) const;

  inline const std::vector<G4float>* ShellProbabilty(std::size_t idx) const;

  void StreamInfo(std::ostream& os) const;

  G4NucLevel(const G4NucLevel &right) = delete;
  G4bool operator==(const G4NucLevel &right) const = delete;
  G4bool operator!=(const G4NucLevel &right) const = delete;
  G4bool operator<(const G4NucLevel &right) const = delete;
  const G4NucLevel& operator=(const G4NucLevel &right) = delete;

private:  

  std::size_t   length;
  G4double fTimeGamma;
  
  std::vector<G4int>    fTrans;
  std::vector<G4float>  fGammaCumProbability;
  std::vector<G4float>  fGammaProbability;
  std::vector<G4float>  fMpRatio;
  std::vector<const std::vector<G4float>*> fShellProbability;
};

inline std::size_t G4NucLevel::NumberOfTransitions() const
{
  return length;
}

inline std::size_t G4NucLevel::FinalExcitationIndex(const std::size_t idx) const
{
  return (std::size_t)(fTrans[idx]/10000);
}

inline G4int G4NucLevel::TransitionType(const std::size_t idx) const
{
  return fTrans[idx]%10000;
}

inline G4double G4NucLevel::GetTimeGamma() const
{
  return fTimeGamma;
}

inline G4float G4NucLevel::GammaProbability(const std::size_t idx) const
{
  return fGammaProbability[idx];
}

inline G4float G4NucLevel::GammaCumProbability(const std::size_t idx) const
{
  return fGammaCumProbability[idx];
}

inline G4float G4NucLevel::MultipolarityRatio(const std::size_t idx) const
{
  return fMpRatio[idx];
}

inline std::size_t G4NucLevel::SampleGammaTransition(const G4double rndm) const
{
  G4float x = rndm;
  std::size_t idx = 0;
  for(; idx<length; ++idx) { 
    if(x <= fGammaCumProbability[idx]) { break; } 
  }
  return idx;
}

inline G4int 
G4NucLevel::SampleShell(const std::size_t idx, const G4double rndm) const
{
  const std::vector<G4float>* prob = fShellProbability[idx];
  G4int i(-1);
  if(nullptr != prob) {
    G4int nn = (G4int)prob->size();
    G4float x = rndm;
    for(i=0; i<nn; ++i) { if(x <= (*prob)[i]) { break; } }
  } 
  return i;
}

inline const std::vector<G4float>* 
G4NucLevel::ShellProbabilty(std::size_t idx) const
{
  return fShellProbability[idx];
}

#endif

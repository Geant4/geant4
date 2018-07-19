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
// $Id: G4UnstableFragmentBreakUp.hh 99812 2016-10-06 08:49:27Z gcosmo $
//
// -------------------------------------------------------------------
//
//      GEANT 4 header file
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4UnstableFragmentBreakUp
//
//      Author:        Vladimir Ivanchenko
//
//      Creation date: 7 May 2010
//
//  Modifications:
// 
// -------------------------------------------------------------------
//  This class providing decay of any fragment on light nucleons using 
//  taking into account only binding energy, for example, it may decay
//  2n -> n + n or 2p -> p + p      
//

#ifndef G4UnstableFragmentBreakUp_h
#define G4UnstableFragmentBreakUp_h 1

#include "globals.hh"
#include "G4VEvaporationChannel.hh"

class G4Fragment;
class G4NuclearLevelData;

class G4UnstableFragmentBreakUp : public G4VEvaporationChannel 
{

public:

  explicit G4UnstableFragmentBreakUp();

  virtual ~G4UnstableFragmentBreakUp();

  virtual G4bool BreakUpChain(G4FragmentVector*, G4Fragment*) final;

  virtual G4double GetEmissionProbability(G4Fragment* fragment) final;

private:

  G4UnstableFragmentBreakUp(const G4UnstableFragmentBreakUp & right) = delete;
  const G4UnstableFragmentBreakUp & operator = 
  (const G4UnstableFragmentBreakUp & right) = delete;
  G4bool operator == (const G4UnstableFragmentBreakUp & right) const = delete;
  G4bool operator != (const G4UnstableFragmentBreakUp & right) const = delete;

  static const G4int Zfr[6];
  static const G4int Afr[6];
  G4double masses[6];

  G4NuclearLevelData* fLevelData;
};

#endif

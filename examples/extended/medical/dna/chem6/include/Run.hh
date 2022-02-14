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
// This example is provided by the Geant4-DNA collaboration
// chem6 example is derived from chem4 and chem5 examples
//
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// J. Appl. Phys. 125 (2019) 104301
// Med. Phys. 45 (2018) e722-e739
// J. Comput. Phys. 274 (2014) 841-882
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157-178
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: W. G. Shin and S. Incerti (CENBG, France)
//
// $Id$
//
/// \file Run.hh
/// \brief Definition of the Run class

#ifndef CHEM6_Run_h
#define CHEM6_Run_h 1

#include "G4Run.hh"
#include "G4THitsMap.hh"

#include "ScoreSpecies.hh"

/// Run class
///
/// In RecordEvent() there is collected information event per event
/// from Hits Collections, and accumulated statistic for the run
class G4VPrimitiveScorer;
class Run : public G4Run
{
public:
  Run();
  virtual ~Run();
  
  virtual void RecordEvent(const G4Event*);
  virtual void Merge(const G4Run*);
  
  G4double GetSumDose()   const { return fSumEne; }
  G4VPrimitiveScorer* GetPrimitiveScorer() const { return fScorerRun;}
  G4THitsMap<G4double>* GetLET() {return fTotalLET;}
  
private:
  G4double fSumEne;
  G4VPrimitiveScorer* fScorerRun;
  G4VPrimitiveScorer* fLETScorerRun;
  G4THitsMap<G4double>* fTotalLET;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

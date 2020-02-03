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
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file Run.hh
/// \brief Definition of the Run class

#ifndef CHEM4_Run_h
#define CHEM4_Run_h 1

#include "ScoreSpecies.hh"
#include "DetectorConstruction.hh"

#include "G4Run.hh"
#include "globals.hh"

/// Run class
///
/// In RecordEvent() there is collected information event per event
/// from Hits Collections, and accumulated statistic for the run
class G4VPrimitiveScorer;
class DetectorConstruction;
class Run : public G4Run
{
public:
  Run();
  virtual ~Run();
  
  virtual void RecordEvent(const G4Event*);
  virtual void Merge(const G4Run*);
  
  G4double GetSumDose()   const { return fSumEne; }
  G4VPrimitiveScorer* GetPrimitiveScorer() const { return fScorerRun;}
  
private:
  G4double fSumEne;
  G4VPrimitiveScorer* fScorerRun;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

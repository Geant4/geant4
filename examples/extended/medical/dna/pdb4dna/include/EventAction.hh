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
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file EventAction.hh
/// \brief Definition of the EventAction class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"

#include <map>

class EventActionMessenger;

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction();

public:
  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

  void AddEdepEvent(G4double edep)
  {
    fTotalEnergyDeposit += edep;
  };
  G4double GetEdepEvent()
  {
    return fTotalEnergyDeposit;
  };

  void AddEdepToNucleotide(G4int numStrand,G4int numNucl,G4double edep)
  {
    if(numStrand==1)
    {
      fEdepStrand1[numNucl]+=edep;
    }
    else{
      fEdepStrand2[numNucl]+=edep;
    }
  }

  void SetEnergyThresForSSB(G4double val)
  {
    fThresEdepForSSB=val;
  };
  void SetDistanceThresForDSB(G4int val)
  {fThresDistForDSB=val;
  };

private:
  // total energy deposit per event
  G4double fTotalEnergyDeposit;
  // map: first strand (G4int : nucleotide ID, G4double : energy deposit)
  std::map<G4int,G4double>  fEdepStrand1;
  // map: second strand (G4int : nucleotide ID, G4double : energy deposit)
  std::map<G4int,G4double>  fEdepStrand2;
  // Min energy to consider single strand break
  G4double fThresEdepForSSB;
  // Max distance to consider double strand break
  G4int fThresDistForDSB;

  EventActionMessenger*     fpEventMessenger;

  // Compute Strand breaks from energy deposits in DNA strands
  void ComputeStrandBreaks(G4int*);
};

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

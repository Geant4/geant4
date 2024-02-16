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
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <iostream>
#include "DetectorConstruction.hh"
#include "G4ThreeVector.hh"

class PrimaryGeneratorAction;
class Run;
class G4Run;
class TrackingAction;

class RunAction : public G4UserRunAction
{
public:
  
  RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
  ~RunAction() override = default;

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

  void  SetRndmFreq(G4int val) { fSaveRndm = val; }
  G4int GetRndmFreq() const { return fSaveRndm; }
  G4Run* GenerateRun() override;   
  
  // Edep in all volume
  G4double GetEdepALL() const { return fEdepAll; }
  void SetEdepALL(G4double vall) { fEdepAll = vall; }
  void AddEdepALL (G4double vall)
  { 
    fEdepAll += vall;
    fEdepAll_err += vall*vall;
  }
  // 0. Edep in homogeneous Medium
  G4double GetEdepMedium() const { return fEdepMedium; }
  void SetEdepMedium(G4double vall) { fEdepMedium = vall; }
  void AddEdepMedium (G4double vall)
  {
    fEdepMedium += vall; 
    fEdepMedium_err += vall*vall;
  }  
  // 1. Edep in Bounding Slice Volume
  G4double GetEdepSlice() const { return fEdepSlice; }
  void SetEdepSlice(G4double vall) { fEdepSlice = vall; }
  void AddEdepSlice (G4double vall)
  { 
    fEdepSlice += vall; 
    fEdepSlice_err += vall*vall;
  }
  // 2. Edep in Soma volume
  G4double GetEdepSoma() const { return fEdepSoma; }
  void SetEdepSoma(G4double vall) { fEdepSoma = vall; }
  void AddEdepSoma (G4double vall)
  { 
   fEdepSoma += vall; 
   fEdepSoma_err += vall*vall;
  }
    
  // 3. Edep in Dendrites volume
  G4double GetEdepDend() const { return fEdepDend; }
  void SetEdepDend(G4double vall){ fEdepDend = vall; }
  void AddEdepDend (G4double vall)
  { 
    fEdepDend += vall; 
    fEdepDend_err += vall*vall;
  }
  
  // 4. Edep in Axon volume
  G4double GetEdepAxon() const { return fEdepAxon; }
  void SetEdepAxon(G4double vall) { fEdepAxon = vall; }
  void AddEdepAxon (G4double vall)
  { 
    fEdepAxon += vall;
    fEdepAxon_err += vall*vall;
  }
  
  // 5. Edep in whole Neuron volume
  G4double GetEdepNeuron() const { return fEdepNeuron; }
  void SetEdepNeuron(G4double vall) { fEdepNeuron = vall; }
  void AddEdepNeuron (G4double vall)
  { 
    fEdepNeuron += vall; 
    fEdepNeuron_err += vall*vall;
  }  
  
  G4int GetNumEvent() const { return fNumEvent; }
  void SetNumEvent(G4int i) { fNumEvent = i; }  

private:

  /////////////////
  // Histogramming
  //
  void CreateHistogram();
  void WriteHistogram();

  /////////////////
  // Print Info
  //
  void PrintRunInfo(const G4Run* run);

  G4bool fDebug{false};

  DetectorConstruction* fDetector;
  PrimaryGeneratorAction* fPrimary;
  Run* fRun{nullptr}; 
  
  //
  // phys 
  G4double fEdepAll{0.0};
  G4double fEdepAll_err{0.0};
  G4double fEdepMedium{0.0};
  G4double fEdepMedium_err{0.0};
  G4double fEdepSlice{0.0};
  G4double fEdepSlice_err{0.0};
  G4double fEdepSoma{0.0};
  G4double fEdepSoma_err{0.0};
  G4double fEdepDend{0.0};
  G4double fEdepDend_err{0.0};
  G4double fEdepAxon{0.0};
  G4double fEdepAxon_err{0.0};
  G4double fEdepNeuron{0.0};
  G4double fEdepNeuron_err{0.0};
  G4int fNumEvent{0};
  G4int fSaveRndm{0};
};

#endif

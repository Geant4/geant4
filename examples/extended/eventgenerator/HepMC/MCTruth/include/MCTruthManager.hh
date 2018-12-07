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
/// \file eventgenerator/HepMC/MCTruth/include/MCTruthManager.hh
/// \brief Definition of the MCTruthManager class
//
//
//
//
// --------------------------------------------------------------
//      GEANT 4 - MCTruthManager class
// --------------------------------------------------------------
//
// Author: Witold POKORSKI (Witold.Pokorski@cern.ch)
// Date  : 2006-02-28
//
// --------------------------------------------------------------
#ifndef INCLUDE_MCTRUTHMANAGER_H 
#define INCLUDE_MCTRUTHMANAGER_H 1

#include "G4Types.hh"
#include "G4LorentzVector.hh"

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

#include "MCTruthConfig.hh"

class MCTruthManager
{

public:

  static MCTruthManager* GetInstance();

  void NewEvent();
  HepMC::GenEvent* GetCurrentEvent() const {return fEvent;}
  void PrintEvent();

  void AddParticle(G4LorentzVector&, G4LorentzVector&, G4LorentzVector&, 
                   G4int, G4int, G4int, G4bool);

  void SetConfig(MCTruthConfig* c) {fConfig=c;}
  MCTruthConfig* GetConfig() const {return fConfig;}

protected:

  MCTruthManager( ); 

  virtual ~MCTruthManager( ); 

private:

  HepMC::GenEvent* fEvent;

  // vector containing barcodes of primary particles (not having any mother)
  //
  std::vector<G4int> fPrimarybarcodes;

  // map containing number of 'segmentations' for each particle (i.e. number
  // of additional vertices introduced in order to attach secondary particles
  // which were created 'in-flight', for instance bremstrahlung gammas, etc)
  //
  std::map<G4int,G4int> fSegmentations;

  // different criteria for storing (or not) particles
  //
  MCTruthConfig* fConfig;

  // recursive printing of the tree
  //
  void PrintTree(HepMC::GenParticle*, G4String);
  
};
#endif // INCLUDE_MCTRUTHMANAGER_H

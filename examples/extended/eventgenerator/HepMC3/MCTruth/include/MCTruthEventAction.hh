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
/// \file eventgenerator/HepMC/MCTruth/include/MCTruthEventAction.hh
/// \brief Definition of the MCTruthEventAction class
//
//
//
//
// --------------------------------------------------------------
//      GEANT 4 - MCTruthEventAction class
// --------------------------------------------------------------
//
// Author: Witold POKORSKI (Witold.Pokorski@cern.ch)
// Date  : 2005-08-29
//
// --------------------------------------------------------------
#ifndef INCLUDE_MCTRUTHEVENTACTION_H 
#define INCLUDE_MCTRUTHEVENTACTION_H 1

#include "G4UserEventAction.hh"
#include "MCTruthManager.hh"

class MCTruthEventAction : public G4UserEventAction
{
public: 

  MCTruthEventAction( ); 

  virtual ~MCTruthEventAction( ); 

  virtual void BeginOfEventAction(const G4Event* anEvent);
  virtual void EndOfEventAction(const G4Event* anEvent);

};

#endif // INCLUDE_MCTRUTHEVENTACTION_H

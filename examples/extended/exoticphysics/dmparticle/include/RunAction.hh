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
/// \file exoticphysics/dmparticle/include/RunAction.hh
/// \brief Definition of the RunAction class
//
// $Id: RunAction.hh 91581 2015-07-27 12:55:57Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   RunAction
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------
//
#ifndef RunAction_h
#define RunAction_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4UserRunAction.hh"
#include "globals.hh"

#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
public: 

  RunAction();
  virtual ~RunAction();

  virtual G4Run* GenerateRun(); 

  virtual void BeginOfRunAction(const G4Run*);
  // In this method histogramms are booked

  virtual void EndOfRunAction(const G4Run*);
  // In this method bookHisto method is called in which histogramms are filled

private:

  // Book predefined histogramms  
  void Book();

private:

  G4AnalysisManager* fAnalysisManager;
  Run*  fRun;
};

#endif


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
// TETRunAction.hh
//
// Author: Haegin Han
// Reference: ICRP Publication 145. Ann. ICRP 49(3), 2020.
// Geant4 Contributors: J. Allison and S. Guatelli
//

#ifndef TETRunAction_h
#define TETRunAction_h 1

#include <ostream>
#include <fstream>
#include <map>

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4UserRunAction.hh"
#include "G4SystemOfUnits.hh"

#include "TETRun.hh"
#include "TETPrimaryGeneratorAction.hh"
#include "TETModelImport.hh"

// *********************************************************************
// The main function of this UserRunAction class is to produce the result
// data and print them.
// -- GenerateRun: Generate TETRun class which will calculate the sum of
//                  energy deposition.
// -- BeginOfRunAction: Set the RunManager to print the progress at the
//                      interval of 10%.
// -- EndOfRunAction: Print the run result by G4cout and std::ofstream.
//  └-- PrintResult: Method to print the result.
// *********************************************************************

class TETRunAction : public G4UserRunAction
{
public:
 explicit TETRunAction(TETModelImport* tetData, G4String output);
 virtual ~TETRunAction() = default;

public:
 virtual G4Run* GenerateRun() override;
 virtual void BeginOfRunAction(const G4Run*) override;
 virtual void EndOfRunAction(const G4Run*) override;

 void PrintResult(std::ostream &out);
  
private:
 TETModelImport* fTetData;
 TETRun*         fRun;
 G4int           fNumOfEvent;
 G4int           fRunID;
 G4String        fOutputFile;
};
#endif






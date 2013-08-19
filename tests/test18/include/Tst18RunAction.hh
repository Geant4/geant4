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
// $Id$
//
//  File:        Tst18RunAction.hh
//  Description: Run action for radioactive decay system test 
//  Author:      F. Lei (DERA UK)
//

#ifndef Tst18RunAction_h
#define Tst18RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>

class Tst18RunActionMessenger;
class G4Run;

class Tst18RunAction : public G4UserRunAction
{
  public:
    Tst18RunAction();
   ~Tst18RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
  
    void SetFilename(G4String val) {fileName = val;};

    void FillParticleName(G4String);
    void FillEnergy(G4double);
    void FillWeight(G4double);
    void FillTime(G4double);
    void FillEventNumber();

  private:
    G4String fileName;  // log file name for the run
    Tst18RunActionMessenger* runMessenger;

    std::vector<G4String> Particles;
    std::vector<G4double> Energies;
    std::vector<G4double> Weights;
    std::vector<G4double> Times;

    G4int EventNumber;

};

#endif


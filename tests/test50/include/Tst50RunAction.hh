//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: Tst50RunAction.hh,v 1.17 2003-07-03 13:43:10 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//
//
// $Id: Tst50RunAction.hh,v 1.17 2003-07-03 13:43:10 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------

#ifndef Tst50RunAction_h
#define Tst50RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class Tst50RunMessenger;
class Tst50AnalysisManager;
class Tst50PrimaryGeneratorAction;
class Tst50DetectorConstruction;
class Tst50RunAction : public G4UserRunAction
{
  public:
    Tst50RunAction();
   ~Tst50RunAction();

  public:
  void BeginOfRunAction(const G4Run*);
  G4int GetRunID ();
  G4bool GetFlag();//returns the choice of test: transmission or 
                   //SP and CSDA range test 
  void EndOfRunAction(const G4Run*);
  void SetTransmissionTest(G4String);
  void EnergyDepositTest(G4String);
  void TransmittedGammaNumber();
  void TransmittedParticleNumber();
  void BackscatteredParticleNumber();

private: 
  Tst50RunMessenger* messenger;
  G4bool flag; // if true transmission test, if false SP and CSDA range
               // valid for massive particles

  G4int runID;
  G4double gammaTransmitted; //number of transmitted gamma 
  G4int numberEvents;// number of events in the BeamOn
  G4double particleTransmitted;// number of transmitted massive particles  
  G4double particleBackscattered;// number of backscattered massive particles 
  Tst50PrimaryGeneratorAction* primary;
  Tst50DetectorConstruction* detector;
  
};
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

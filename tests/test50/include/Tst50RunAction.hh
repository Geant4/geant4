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
// $Id: Tst50RunAction.hh,v 1.15 2003-05-17 14:18:48 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//
//
// $Id: Tst50RunAction.hh,v 1.15 2003-05-17 14:18:48 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 02 Mar  1999   S. Guatelli   1st implementation
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
  void Set_Trans(G4String);
  G4bool Get_flag();
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
  void Trans_number();
  void Back_number();
  void gamma_transmitted();
  G4int GetRun_ID ();

private:
  G4double particle_trans;// number of transmitted massive particles
  G4double particle_back;// number of backscattered massive particles
  G4int numberEvents;// number of events in the BeamOn 
  G4int runID;
  G4double gamma_trans;// number of transmitted gamma 
  Tst50PrimaryGeneratorAction* p_Primary;
  Tst50DetectorConstruction* p_Detector;
  G4bool flag; // if true transmission test, if false SP and CSDA range
               // valid for massive particles
  Tst50RunMessenger* p_messenger;
};
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
//
// $Id: Tst50SteppingAction.hh,v 1.16 2006-06-29 22:05:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 27 May  2003   S.Guatelli    first code review 
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------

#ifndef Tst50SteppingAction_h
#define Tst50SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

class G4Step;
class Tst50AnalysisManager;
class Tst50RunAction;
class Tst50PrimaryGeneratorAction;
class Tst50DetectorConstruction;
class Tst50SteppingAction : public G4UserSteppingAction
{
public:

  Tst50SteppingAction(Tst50PrimaryGeneratorAction*,
		      Tst50RunAction*, 
		      Tst50DetectorConstruction*);

  ~Tst50SteppingAction();

  void UserSteppingAction(const G4Step* aStep);

private:

  Tst50PrimaryGeneratorAction* primaryAction;
  Tst50RunAction* runAction; 
  Tst50DetectorConstruction* detector;     
};
#endif

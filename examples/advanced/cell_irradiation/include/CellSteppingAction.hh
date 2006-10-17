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
//    **************************************
//    *                                    *
//    *        CellSteppingAction.hh       *
//    *                                    *
//    **************************************
//
// 
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//	   Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
//
// History:
// -----------
// 20 September 2006   S. Guatelli, B. Mascialino   1st implementation
//
// -------------------------------------------------------------------

#ifndef CellSteppingAction_h
#define CellSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

class G4Step;
class CellAnalysisManager;
class CellRunAction;
class CellPrimaryGeneratorAction;
class CellDetectorConstruction;
class CellSteppingAction : public G4UserSteppingAction
{
public:

   CellSteppingAction(CellPrimaryGeneratorAction*,
		      CellRunAction*, 
		      CellDetectorConstruction*);

  ~CellSteppingAction();

  void UserSteppingAction(const G4Step* aStep);

private:

  CellPrimaryGeneratorAction* primaryAction;
  CellRunAction* runAction; 
  CellDetectorConstruction* detector;     
};
#endif

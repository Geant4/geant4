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
#ifndef PARN02JOB_HH
#define PARN02JOB_HH
#include "G4VtbbJob.hh"

class G4tbbRunManager;
class ExN02DetectorConstruction;

class ParN02Job : public G4VtbbJob {
public:
  ParN02Job(const G4String& macroFile);
  ~ParN02Job();
  //Interface implementation, from base class
  void CreateDetector(G4tbbRunManager* rm);
  void UserActions(G4tbbRunManager* rm);
  void InitSetup(G4tbbRunManager* rm );
  void JobPrepare(G4tbbRunManager* rm );
private:
  ExN02DetectorConstruction* detector;
};

#endif

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
// $Id: TiaraSim.hh,v 1.5 2006/06/29 15:44:14 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// ----------------------------------------------------------------------
//
// Class TiaraSim
//

#ifndef TiaraSim_hh
#define TiaraSim_hh TiaraSim_hh
#include "G4UIterminal.hh"
#include <map>
#include <string>
#include "G4ProcessPlacer.hh"

class G4RunManager;
class G4VUserPrimaryGeneratorAction;
class G4VUserPhysicsList;
class TiaraPhysicsList;
class G4VUserDetectorConstruction;
class G4VisManager;
class TiaraStackingAction;
class G4UserEventAction;
class TiaraEnergyCutProcess;

class TiaraSim {
public:
  ~TiaraSim();
  static TiaraSim &GetTiaraSim();
  void SetGeometry(G4VUserDetectorConstruction *geo);
  void SetPrimaryGenerator(G4VUserPrimaryGeneratorAction *primGen);
  void SetPhysicsList(G4VUserPhysicsList *physList);
  void initialize();
  void initializeVisManager();
  void startSession();
  void AddParticleCut(const std::string &particle,
		      G4double cut);
  void AddTiaraEventAction(G4UserEventAction *eventAction);
  void AddVisRunAction();
  void BeamOn(G4int nEvents);

private:
  TiaraSim();
  static TiaraSim *fTiaraSim;
  G4RunManager *frunMgr;
  G4VisManager *fVisManager;
  G4UIterminal fSession;

  G4VUserDetectorConstruction *fGeometry;
  G4VUserPrimaryGeneratorAction *fPrimary;
  G4VUserPhysicsList *fPhysics;

  std::vector<TiaraEnergyCutProcess *> fCutProcessVector;
  std::vector<G4ProcessPlacer> fPlacers;

  TiaraStackingAction *fStackingAction;
  G4UserEventAction *fUserEventAction;
};




#endif











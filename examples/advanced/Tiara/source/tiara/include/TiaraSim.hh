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
// $Id: TiaraSim.hh,v 1.4 2003/06/25 09:12:51 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
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











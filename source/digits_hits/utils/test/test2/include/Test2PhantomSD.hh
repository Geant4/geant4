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

#ifndef Test2PhantomSD_h
#define Test2PhantomSD_h 1

#include "G4VSensitiveDetector.hh"
#include "Test2PhantomHit.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"
#include "G4GeometryTolerance.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class G4VSolid;

class Test2PhantomSD : public G4VSensitiveDetector {

public:
  Test2PhantomSD(G4String name, G4int segment[3]);
  ~Test2PhantomSD();

  void Initialize(G4HCofThisEvent * HCE);
  G4bool ProcessHits(G4Step * aStep, G4TouchableHistory * ROhist);
  void EndOfEvent(G4HCofThisEvent * HCE);
  void clear();
  void DrawAll();
  void PrintAll();

private:
  G4VSolid* GetSolid(G4Step* aStep);
  G4double GetVolume(G4Step* aStep);
  G4double GetArea(G4Step* aStep);
  G4int IsSelectedSurface(G4Step* aStep);
  G4bool IsPassed(G4Step* aStep);
  G4double GetAngleFactor(G4Step* aStep,G4int dirFlag);
  G4bool IsSecondary(G4Step* aStep);
  G4bool IsEnterOrFirstStep(G4Step* aStep);  
  G4bool IsExit(G4Step* aStep);

private:
  Test2PhantomHitsCollection * fPhantomCollection;
  G4int nSegment[3];

  G4int fCurrentTrkID;
  G4double fCellTrack;

};

#endif


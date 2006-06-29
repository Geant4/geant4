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

#ifndef PhotInStackingAction_H
#define PhotInStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

#include "PhotInConstants.hh"

#include "G4Track.hh"
//#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

class PhotInStackingAction : public G4UserStackingAction
{
public:
  PhotInStackingAction();
  virtual ~PhotInStackingAction();

  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
  virtual void PrepareNewEvent();

  static G4int GetNGamma     (G4int i) { return PhotInNGamma     [i]; }
  static G4int GetNElectron  (G4int i) { return PhotInNElectron  [i]; }
  static G4int GetNPositron  (G4int i) { return PhotInNPositron  [i]; }
  static G4int GetNNeutrons  (G4int i) { return PhotInNNeutrons  [i]; }
  static G4int GetNProtons   (G4int i) { return PhotInNProtons   [i]; }
  static G4int GetNDeuterons (G4int i) { return PhotInNDeuterons [i]; }
  static G4int GetNTritons   (G4int i) { return PhotInNTritons   [i]; }
  static G4int GetNHe3s    		(G4int i) { return PhotInNHe3s    		[i]; }
  static G4int GetNAlphas    (G4int i) { return PhotInNAlphas    [i]; }
  static G4int GetNLambdas   (G4int i) { return PhotInNLambdas   [i]; }
  static G4int GetNHeavyFrags(G4int i) { return PhotInNHeavyFrags[i]; }
  static G4int GetNMesons    (G4int i) { return PhotInNMesons    [i]; }

  static G4double GetEMinGamma     (G4int i) { return PhotInEMinGamma     [i]; }
  static G4double GetEMinElectron  (G4int i) { return PhotInEMinElectron  [i]; }
  static G4double GetEMinPositron  (G4int i) { return PhotInEMinPositron  [i]; }
  static G4double GetEMinNeutrons  (G4int i) { return PhotInEMinNeutrons  [i]; }
  static G4double GetEMinProtons   (G4int i) { return PhotInEMinProtons   [i]; }
  static G4double GetEMinDeuterons (G4int i) { return PhotInEMinDeuterons [i]; }
  static G4double GetEMinTritons   (G4int i) { return PhotInEMinTritons   [i]; }
  static G4double GetEMinHe3s      (G4int i) { return PhotInEMinHe3s      [i]; }
  static G4double GetEMinAlphas    (G4int i) { return PhotInEMinAlphas    [i]; }
  static G4double GetEMinLambdas   (G4int i) { return PhotInEMinLambdas   [i]; }
  static G4double GetEMinHeavyFrags(G4int i) { return PhotInEMinHeavyFrags[i]; }
  static G4double GetEMinMesons    (G4int i) { return PhotInEMinMesons    [i]; }

  static G4double GetEMaxGamma     (G4int i) { return PhotInEMaxGamma     [i]; }
  static G4double GetEMaxElectron  (G4int i) { return PhotInEMaxElectron  [i]; }
  static G4double GetEMaxPositron  (G4int i) { return PhotInEMaxPositron  [i]; }
  static G4double GetEMaxNeutrons  (G4int i) { return PhotInEMaxNeutrons  [i]; }
  static G4double GetEMaxProtons   (G4int i) { return PhotInEMaxProtons   [i]; }
  static G4double GetEMaxDeuterons (G4int i) { return PhotInEMaxDeuterons [i]; }
  static G4double GetEMaxTritons   (G4int i) { return PhotInEMaxTritons   [i]; }
  static G4double GetEMaxHe3s      (G4int i) { return PhotInEMaxHe3s      [i]; }
  static G4double GetEMaxAlphas    (G4int i) { return PhotInEMaxAlphas    [i]; }
  static G4double GetEMaxLambdas   (G4int i) { return PhotInEMaxLambdas   [i]; }
  static G4double GetEMaxHeavyFrags(G4int i) { return PhotInEMaxHeavyFrags[i]; }
  static G4double GetEMaxMesons    (G4int i) { return PhotInEMaxMesons    [i]; }

private: //--- BODY ---
  // Can not be moved to the PhotInConstants.hh as static because there is a vector of Hits
  static G4int PhotInNGamma     [PhotInDiNSections];
  static G4int PhotInNElectron  [PhotInDiNSections];
  static G4int PhotInNPositron  [PhotInDiNSections];
  static G4int PhotInNNeutrons  [PhotInDiNSections];
  static G4int PhotInNProtons   [PhotInDiNSections];
  static G4int PhotInNDeuterons [PhotInDiNSections];
  static G4int PhotInNTritons   [PhotInDiNSections];
  static G4int PhotInNHe3s      [PhotInDiNSections];		
  static G4int PhotInNAlphas    [PhotInDiNSections];
  static G4int PhotInNLambdas   [PhotInDiNSections];
  static G4int PhotInNHeavyFrags[PhotInDiNSections];//Includes any nuclearFragment with A>4
  static G4int PhotInNMesons    [PhotInDiNSections];//Includes pions and kaons of all signs

  static G4double PhotInEMinGamma     [PhotInDiNSections];
  static G4double PhotInEMinElectron  [PhotInDiNSections];
  static G4double PhotInEMinPositron  [PhotInDiNSections];
  static G4double PhotInEMinNeutrons  [PhotInDiNSections];
  static G4double PhotInEMinProtons   [PhotInDiNSections];
  static G4double PhotInEMinDeuterons [PhotInDiNSections];
  static G4double PhotInEMinTritons   [PhotInDiNSections];
  static G4double PhotInEMinHe3s      [PhotInDiNSections];
  static G4double PhotInEMinAlphas    [PhotInDiNSections];
  static G4double PhotInEMinLambdas   [PhotInDiNSections];
  static G4double PhotInEMinHeavyFrags[PhotInDiNSections];
  static G4double PhotInEMinMesons    [PhotInDiNSections];

  static G4double PhotInEMaxGamma     [PhotInDiNSections];
  static G4double PhotInEMaxElectron  [PhotInDiNSections];
  static G4double PhotInEMaxPositron  [PhotInDiNSections];
  static G4double PhotInEMaxNeutrons  [PhotInDiNSections];
  static G4double PhotInEMaxProtons   [PhotInDiNSections];
  static G4double PhotInEMaxDeuterons [PhotInDiNSections];
  static G4double PhotInEMaxTritons   [PhotInDiNSections];
  static G4double PhotInEMaxHe3s      [PhotInDiNSections];
  static G4double PhotInEMaxAlphas    [PhotInDiNSections];
  static G4double PhotInEMaxLambdas   [PhotInDiNSections];
  static G4double PhotInEMaxHeavyFrags[PhotInDiNSections];
  static G4double PhotInEMaxMesons    [PhotInDiNSections];
};									

#endif


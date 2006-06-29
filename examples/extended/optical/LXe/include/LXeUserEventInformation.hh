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
#include "G4VUserEventInformation.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#ifndef LXeUserEventInformation_h
#define LXeUserEventInformation_h 1

class LXeUserEventInformation : public G4VUserEventInformation
{
public:
  LXeUserEventInformation();
  ~LXeUserEventInformation();
  
  inline void Print()const{};

  void IncPhotonCount_Scint(){photonCount_Scint++;}
  void IncPhotonCount_Ceren(){photonCount_Ceren++;}
  void IncEDep(G4double dep){totE+=dep;}
  void IncAbsorption(){absorptionCount++;}
  void IncBoundaryAbsorption(){boundaryAbsorptionCount++;}
  void IncHitCount(G4int i=1){hitCount+=i;}

  void SetEWeightPos(const G4ThreeVector& p){eWeightPos=p;}
  void SetReconPos(const G4ThreeVector& p){reconPos=p;}
  void SetConvPos(const G4ThreeVector& p){convPos=p;convPosSet=true;}
  void SetPosMax(const G4ThreeVector& p,G4double edep){posMax=p;edepMax=edep;}

  G4int GetPhotonCount_Scint()const {return photonCount_Scint;}
  G4int GetPhotonCount_Ceren()const {return photonCount_Ceren;}
  G4int GetHitCount()const {return hitCount;}
  G4double GetEDep()const {return totE;}
  G4int GetAbsorptionCount()const {return absorptionCount;}
  G4int GetBoundaryAbsorptionCount() const {return boundaryAbsorptionCount;}
  
  G4ThreeVector GetEWeightPos(){return eWeightPos;}
  G4ThreeVector GetReconPos(){return reconPos;}
  G4ThreeVector GetConvPos(){return convPos;}
  G4ThreeVector GetPosMax(){return posMax;}
  G4double GetEDepMax(){return edepMax;}
  G4double IsConvPosSet(){return convPosSet;}

  //Gets the total photon count produced
  G4int GetPhotonCount(){return photonCount_Scint+photonCount_Ceren;}

  void IncPMTSAboveThreshold(){pmtsAboveThreshold++;}
  G4int GetPMTSAboveThreshold(){return pmtsAboveThreshold;}

private:

  G4int hitCount;
  G4int photonCount_Scint;
  G4int photonCount_Ceren;
  G4int absorptionCount;
  G4int boundaryAbsorptionCount;

  G4double totE;  

  //These only have meaning if totE > 0
  //If totE = 0 then these wont be set by EndOfEventAction
  G4ThreeVector eWeightPos;
  G4ThreeVector reconPos; //Also relies on hitCount>0
  G4ThreeVector convPos;//true (initial) converstion position
  G4bool convPosSet;
  G4ThreeVector posMax;
  G4double edepMax;

  G4int pmtsAboveThreshold;

};

#endif






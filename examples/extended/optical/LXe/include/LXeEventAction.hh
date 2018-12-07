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
/// \file optical/LXe/include/LXeEventAction.hh
/// \brief Definition of the LXeEventAction class
//

#ifndef LXeEventAction_h
#define LXeEventAction_h 1

#include "LXeEventMessenger.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4Event;
class LXeDetectorConstruction;

class LXeEventAction : public G4UserEventAction
{
  public:

    LXeEventAction(const LXeDetectorConstruction*);
    virtual ~LXeEventAction();

  public:

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    void SetEventVerbose(G4int v){fVerbose=v;}

    void SetPMTThreshold(G4int t){fPMTThreshold=t;}

    void SetForceDrawPhotons(G4bool b){fForcedrawphotons=b;}
    void SetForceDrawNoPhotons(G4bool b){fForcenophotons=b;}

    void IncPhotonCount_Scint(){fPhotonCount_Scint++;}
    void IncPhotonCount_Ceren(){fPhotonCount_Ceren++;}
    void IncEDep(G4double dep){fTotE+=dep;}
    void IncAbsorption(){fAbsorptionCount++;}
    void IncBoundaryAbsorption(){fBoundaryAbsorptionCount++;}
    void IncHitCount(G4int i=1){fHitCount+=i;}

    void SetEWeightPos(const G4ThreeVector& p){fEWeightPos=p;}
    void SetReconPos(const G4ThreeVector& p){fReconPos=p;}
    void SetConvPos(const G4ThreeVector& p){fConvPos=p;fConvPosSet=true;}
    void SetPosMax(const G4ThreeVector& p,G4double edep) {
      fPosMax = p;
      fEdepMax = edep;
    }

    G4int GetPhotonCount_Scint()const {return fPhotonCount_Scint;}
    G4int GetPhotonCount_Ceren()const {return fPhotonCount_Ceren;}
    G4int GetHitCount()const {return fHitCount;}
    G4double GetEDep()const {return fTotE;}
    G4int GetAbsorptionCount()const {return fAbsorptionCount;}
    G4int GetBoundaryAbsorptionCount() const {return fBoundaryAbsorptionCount;}

    G4ThreeVector GetEWeightPos(){return fEWeightPos;}
    G4ThreeVector GetReconPos(){return fReconPos;}
    G4ThreeVector GetConvPos(){return fConvPos;}
    G4ThreeVector GetPosMax(){return fPosMax;}
    G4double GetEDepMax(){return fEdepMax;}
    G4double IsConvPosSet(){return fConvPosSet;}

    //Gets the total photon count produced
    G4int GetPhotonCount(){return fPhotonCount_Scint+fPhotonCount_Ceren;}

    void IncPMTSAboveThreshold(){fPMTsAboveThreshold++;}
    G4int GetPMTSAboveThreshold(){return fPMTsAboveThreshold;}

  private:

    LXeEventMessenger* fEventMessenger;
    const LXeDetectorConstruction* fDetector;

    G4int fScintCollID;
    G4int fPMTCollID;

    G4int fVerbose;

    G4int fPMTThreshold;

    G4bool fForcedrawphotons;
    G4bool fForcenophotons;

    G4int fHitCount;
    G4int fPhotonCount_Scint;
    G4int fPhotonCount_Ceren;
    G4int fAbsorptionCount;
    G4int fBoundaryAbsorptionCount;

    G4double fTotE;

    //These only have meaning if totE > 0
    //If totE = 0 then these wont be set by EndOfEventAction
    G4ThreeVector fEWeightPos;
    G4ThreeVector fReconPos; //Also relies on hitCount>0
    G4ThreeVector fConvPos;//true (initial) converstion position
    G4bool fConvPosSet;
    G4ThreeVector fPosMax;
    G4double fEdepMax;

    G4int fPMTsAboveThreshold;
};

#endif

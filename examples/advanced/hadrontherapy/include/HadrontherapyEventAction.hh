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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#ifndef HadrontherapyEventAction_h
#define HadrontherapyEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class HadrontherapyMatrix;
class HadrontherapyEventActionMessenger;

class HadrontherapyEventAction : public G4UserEventAction
{
public:
  HadrontherapyEventAction();
  ~HadrontherapyEventAction();

public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);

  void SetPrintModulo(G4int val)
  {
    printModulo = val;
  };

  void SetDrawFlag(G4String val)
  {
    drawFlag = val;
  };

private: 
  G4String drawFlag; //Visualisation flag
  G4int hitsCollectionID;
  //HadrontherapyMatrix *matrix; 
  G4int printModulo;  
  HadrontherapyEventActionMessenger* pointerEventMessenger;
};

class B4aEventAction : public G4UserEventAction
{
public:
    B4aEventAction();
    virtual ~B4aEventAction();
    
    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
    
    void AddAbs(G4double de, G4double dl);
    void AddGap(G4double de, G4double dl);
    
private:
    G4double  fEnergyAbs;
    G4double  fEnergyGap;
    G4double  fTrackLAbs;
    G4double  fTrackLGap;
};

// inline functions

inline void B4aEventAction::AddAbs(G4double de, G4double dl) {
    fEnergyAbs += de;
    fTrackLAbs += dl;
}

inline void B4aEventAction::AddGap(G4double de, G4double dl) {
    fEnergyGap += de;
    fTrackLGap += dl;
}


#endif




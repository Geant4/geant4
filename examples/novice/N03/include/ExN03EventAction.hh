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
// $Id: ExN03EventAction.hh,v 1.10 2006/06/29 17:48:36 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExN03EventAction_h
#define ExN03EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class ExN03RunAction;
class ExN03EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN03EventAction : public G4UserEventAction
{
 public:
   ExN03EventAction(ExN03RunAction*);
  ~ExN03EventAction();

 public:
   void  BeginOfEventAction(const G4Event*);
   void    EndOfEventAction(const G4Event*);
    
   void AddAbs(G4double de, G4double dl) {EnergyAbs += de; TrackLAbs += dl;};
   void AddGap(G4double de, G4double dl) {EnergyGap += de; TrackLGap += dl;};
                     
   void SetDrawFlag   (G4String val)  {drawFlag = val;};
   void SetPrintModulo(G4int    val)  {printModulo = val;};
    
 private:
   ExN03RunAction*  runAct;
   
   G4double  EnergyAbs, EnergyGap;
   G4double  TrackLAbs, TrackLGap;
                     
   G4String  drawFlag;
   G4int     printModulo;
                             
   ExN03EventActionMessenger*  eventMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    

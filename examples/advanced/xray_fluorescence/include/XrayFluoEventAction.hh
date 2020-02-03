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
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//
// -------------------------------------------------------------------


#ifndef XrayFluoEventAction_h
#define XrayFluoEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "XrayFluoDetectorConstruction.hh"
#include "XrayFluoPlaneDetectorConstruction.hh"
#include "XrayFluoMercuryDetectorConstruction.hh"

class XrayFluoRunAction;
class XrayFluoEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


class XrayFluoEventAction : public G4UserEventAction
{
public:
  
  XrayFluoEventAction(const XrayFluoDetectorConstruction*);
  XrayFluoEventAction(const XrayFluoPlaneDetectorConstruction*);
  XrayFluoEventAction(const XrayFluoMercuryDetectorConstruction*);

  virtual ~XrayFluoEventAction();
  
public:
  virtual void   BeginOfEventAction(const G4Event*);
  virtual void   EndOfEventAction(const G4Event*);
  
  //method used to set the flags for drawing the tracks
  void SetDrawFlag   (G4String val)  {drawFlag = val;};
  void SetPrintModulo(G4int    val)  {printModulo = val;};
  
private:
  
  G4String                    drawFlag;
  G4int                       HPGeCollID; 
  
  //pointer to XrayFluoEventActionMessenger
  XrayFluoEventActionMessenger*  eventMessenger;
  G4int                       printModulo;                         
  
  //this method generates a gaussian distribution
  //the argument is the mean
  //the sigma is to be set in the file XrayFluoEventAction.cc 
  G4double RandomCut(G4double);

  XrayFluoVDetectorType* detectorType;

  
   //XrayFluoRunAction* runManager;


  //this method distributes the energy deposit (which must be given as
  //argument) according to the response function stored in the file 
  //response.dat
  //public: G4double ResponseFunction(G4double);

  

};

#endif

    

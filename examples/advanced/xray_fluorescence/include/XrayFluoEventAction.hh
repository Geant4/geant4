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
//
// $Id: XrayFluoEventAction.hh
// GEANT4 tag $Name: xray_fluo-V04-01-03
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

class XrayFluoRunAction;
class XrayFluoEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


class XrayFluoEventAction : public G4UserEventAction
{
public:
  
  XrayFluoEventAction();
  
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

  //this method distributes the energy deposit (which must be given as
  //argument) according to the response function stored in the file 
  //response.dat
  
  
  XrayFluoRunAction* runManager;

public: G4double ResponseFunction(G4double);
  

};

#endif

    

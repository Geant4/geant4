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
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//  19 Jun 2003  Alfonso Mantero Created
//
// -------------------------------------------------------------------

#ifndef XrayFluoVDetectorType_hh
#define XrayFluoVDetectorType_hh 1

#include "globals.hh"

class XrayFluoVDetectorType

{
public:
  
  virtual ~XrayFluoVDetectorType();
  virtual G4String GetDetectorMaterial()=0;
  virtual G4double ResponseFunction(G4double)=0;
  //returns a random value of the energy measured according to the tabulated
  //peack just lower of the energy deposited. the third integer not always is used. 
  //when it used, user must specify it, otherwise, he gets a runtime error.
  virtual G4double GetInfData(G4double,G4double,G4int=0)=0;
  //returns a random value of the energy measured according to the tabulated
  //peack just upper of the energy deposited. the third integer not always is used. 
  //when it used, user must specify it, otherwise, he gets a runtime error.
  virtual G4double GetSupData(G4double,G4double,G4int=0)=0;
  virtual void LoadResponseData(G4String)=0;
  virtual void LoadEfficiencyData(G4String)=0;
  
protected:
  
  XrayFluoVDetectorType();
};

#endif

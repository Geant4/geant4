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
// $Id: XrayFluoVdetectorType.hh
// GEANT4 tag $Name:
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

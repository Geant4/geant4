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
// 28 Nov 2001 Elena Guardincerri     Created
// 29 Nov 2002 Alfonso Mantero minor modifications (Alfonso.mantero@ge.infn.it)
// -------------------------------------------------------------------

#include "XrayFluoSensorHit.hh"

G4ThreadLocal G4Allocator<XrayFluoSensorHit> *XrayFluoSensorHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoSensorHit::XrayFluoSensorHit()
{
  EdepTot = 0.;  
  EdepDetect = 0.;

 
  //  G4cout << "XrayFluoSensorHit created" << G4endl;   
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoSensorHit::XrayFluoSensorHit(const XrayFluoSensorHit& right):
  G4VHit(right)
{
  EdepTot = right.EdepTot ; 

  EdepDetect = right.EdepDetect;

  //G4cout << "XrayFluoSensorHit created" << G4endl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoSensorHit::~XrayFluoSensorHit()
{
  //G4cout << "XrayFluo SensorHit deleted" << G4endl; 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const XrayFluoSensorHit& XrayFluoSensorHit::operator=(const XrayFluoSensorHit& right)
{
  EdepTot = right.EdepTot ; 
  EdepDetect = right.EdepDetect;
 
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool XrayFluoSensorHit::operator==(const XrayFluoSensorHit&) const
{
  return false;
}


void XrayFluoSensorHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoSensorHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....











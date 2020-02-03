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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "XrayFluoActionInitializer.hh"
#include "XrayFluoDetectorConstruction.hh"
#include "XrayFluoPlaneDetectorConstruction.hh"
#include "XrayFluoMercuryDetectorConstruction.hh"
#include "XrayFluoPrimaryGeneratorAction.hh"
#include "XrayFluoPlanePrimaryGeneratorAction.hh"
#include "XrayFluoMercuryPrimaryGeneratorAction.hh"
#include "XrayFluoRunAction.hh"
#include "XrayFluoEventAction.hh"
#include "XrayFluoSteppingAction.hh"
#include "XrayFluoSteppingVerbose.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XrayFluoActionInitializer::XrayFluoActionInitializer(G4int geometryFlag) : 
  G4VUserActionInitialization(),fGeometryFlag(geometryFlag)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XrayFluoActionInitializer::Build() const 
{
  //Actual version will depend on the geometry
  XrayFluoEventAction* eventAction = 0;
  XrayFluoSteppingAction* stepAction = new XrayFluoSteppingAction();

   //Selecting the PrimaryGenerator depending upon Geometry setup selected
  if (fGeometryFlag == 1 || fGeometryFlag == 4) 
    {  
      const XrayFluoDetectorConstruction* testBeamDetector = 
	static_cast<const XrayFluoDetectorConstruction*>
	(G4RunManager::GetRunManager()->GetUserDetectorConstruction());  
      
      eventAction = new XrayFluoEventAction(testBeamDetector);
      SetUserAction(new XrayFluoPrimaryGeneratorAction(testBeamDetector));
    }
  else if (fGeometryFlag == 2) 
    {      
      const XrayFluoPlaneDetectorConstruction* planeDetector = 
	static_cast<const XrayFluoPlaneDetectorConstruction*>
	(G4RunManager::GetRunManager()->GetUserDetectorConstruction());  

      eventAction = new XrayFluoEventAction(planeDetector);
      SetUserAction(new XrayFluoPlanePrimaryGeneratorAction(planeDetector));
    }
  else if (fGeometryFlag == 3) 
    {      
      const XrayFluoMercuryDetectorConstruction* mercuryDetector = 
	static_cast<const XrayFluoMercuryDetectorConstruction*>
	(G4RunManager::GetRunManager()->GetUserDetectorConstruction());  
      eventAction = new XrayFluoEventAction(mercuryDetector);
      stepAction->SetMercuryFlag(true);
      SetUserAction(new XrayFluoMercuryPrimaryGeneratorAction(mercuryDetector));
  }
  SetUserAction(eventAction);
  SetUserAction(new XrayFluoRunAction());
  SetUserAction(stepAction);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XrayFluoActionInitializer::BuildForMaster() const
{
  SetUserAction(new XrayFluoRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4VSteppingVerbose* XrayFluoActionInitializer::InitializeSteppingVerbose() const
{
  //XrayFluo Verbose output class
  return new XrayFluoSteppingVerbose();
}

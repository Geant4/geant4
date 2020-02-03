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

#include "TransitionRadiationPhysics.hh"
#include "DetectorConstruction.hh"
#include "G4VXTRenergyLoss.hh"
#include "G4ProcessManager.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4VXTRenergyLoss.hh"
#include "G4RegularXTRadiator.hh"
#include "G4TransparentRegXTRadiator.hh"
#include "G4GammaXTRadiator.hh"
#include "G4StrawTubeXTRadiator.hh"

#include "G4XTRGammaRadModel.hh"
#include "G4XTRRegularRadModel.hh"
#include "G4XTRTransparentRegRadModel.hh"
#include "XTRTransparentRegRadModel.hh"

G4ThreadLocal 
G4VXTRenergyLoss* TransitionRadiationPhysics::fXTRProcess = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TransitionRadiationPhysics::TransitionRadiationPhysics(G4int verb,
    DetectorConstruction* ptr) 
  : G4VPhysicsConstructor("XTR"), 
    fDetector(ptr), 
    fVerbose(verb),
    fXTRModel("transpM")
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TransitionRadiationPhysics::~TransitionRadiationPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TransitionRadiationPhysics::ConstructProcess()
{
  if("dummy" == fXTRModel) { return; }
  if(0 < fVerbose) {
    G4cout<< "TransitionRadiationPhysics: XTR model <" << fXTRModel
          << ">" <<G4endl;
  }
  RadiatorDescription* rDescription = fDetector->GetRadiatorDescription();

  if(fXTRModel == "gammaR" ) {      

    fXTRProcess = new G4GammaXTRadiator(rDescription->fLogicalVolume,
                                        100., 100.,
                                        rDescription->fFoilMaterial,
                                        rDescription->fGasMaterial,
                                        rDescription->fFoilThickness,
                                        rDescription->fGasThickness,
                                        rDescription->fFoilNumber,
                                        "GammaXTRadiator");
  }
  else if(fXTRModel == "gammaM" ) 
  {
    fXTRProcess = new G4XTRGammaRadModel(rDescription->fLogicalVolume,
                                         100., 100.,
                                         rDescription->fFoilMaterial,
                                         rDescription->fGasMaterial,
                                         rDescription->fFoilThickness,
                                         rDescription->fGasThickness,
                                         rDescription->fFoilNumber,
                                         "GammaXTRadiator");
  }
  else if(fXTRModel == "strawR" ) 
  {
    fXTRProcess = new G4StrawTubeXTRadiator(rDescription->fLogicalVolume,
                                            rDescription->fFoilMaterial,
                                            rDescription->fGasMaterial,
                                            0.53,
                                            3.14159,
                                            fDetector->GetAbsorberMaterial(),
                                            true,
                                            "strawXTRadiator");
  }
  else if(fXTRModel == "regR" ) 
  {      
    fXTRProcess = new G4RegularXTRadiator(rDescription->fLogicalVolume,
                                          rDescription->fFoilMaterial,
                                          rDescription->fGasMaterial,
                                          rDescription->fFoilThickness,
                                          rDescription->fGasThickness,
                                          rDescription->fFoilNumber,
                                          "RegularXTRadiator");            
  }
  else if(fXTRModel == "transpR" ) 
  {
    // G4TransparentRegXTRadiator* 
    fXTRProcess = new G4TransparentRegXTRadiator(rDescription->fLogicalVolume,
                                         rDescription->fFoilMaterial,
                                         rDescription->fGasMaterial,
                                         rDescription->fFoilThickness,
                                         rDescription->fGasThickness,
                                         rDescription->fFoilNumber,
                                         "RegularXTRadiator");
  }
  else if(fXTRModel == "regM" ) 
  {
    fXTRProcess = new G4XTRRegularRadModel(rDescription->fLogicalVolume,
                                           rDescription->fFoilMaterial,
                                           rDescription->fGasMaterial,
                                           rDescription->fFoilThickness,
                                           rDescription->fGasThickness,
                                           rDescription->fFoilNumber,
                                           "RegularXTRadiator");
       
  }
  else if(fXTRModel == "transpM" ) 
  { 
    fXTRProcess = new XTRTransparentRegRadModel(rDescription->fLogicalVolume,
                                                rDescription->fFoilMaterial,
                                                rDescription->fGasMaterial,
                                                rDescription->fFoilThickness,
                                                rDescription->fGasThickness,
                                                rDescription->fFoilNumber,
                                                "RegularXTRadiator");
  }     
  if(!fXTRProcess) { 
    if(0 < fVerbose) {
      G4cout<< "TransitionRadiationPhysics: XTR model <" << fXTRModel
            << "> is not known - no XTR process defined" <<G4endl;
    }
    return; 
  }

  fXTRProcess->SetVerboseLevel(fVerbose);

  G4Electron* elec = G4Electron::Electron();
  G4ProcessManager* manager = elec->GetProcessManager();
  manager->AddDiscreteProcess(fXTRProcess);

  G4Positron* posi = G4Positron::Positron();
  manager = posi->GetProcessManager();
  manager->AddDiscreteProcess(fXTRProcess);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


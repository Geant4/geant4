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
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            CompPrimaryGeneratorAction.cc                *
// * -------                                                            *
// *                                                                    *
// * Version:           0.5                                             *
// * Date:              15/11/00                                        *
// * Author:            F.Lei                                           *
// * Organisation:      DERA, UK                                        *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 15.11.2000 F.Lei
// - New version of PrimaryGeneratorAction using the GPS insead of the
//   standard particle gun.  
// - The PrimaryGeneratorMessenger is no longer needed.
// 
// 06.11.2000 R.Nartallo
// - First implementation of PrimaryGeneratorAction
// - Based on Chandra and XMM models by S Magni and F Lei
//
// **********************************************************************

#include "ZIIIPrimaryGeneratorAction.hh"

#include "ZIIIDetectorConstruction.hh"


#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//ZIIIPrimaryGeneratorAction::ZIIIPrimaryGeneratorAction()

ZIIIPrimaryGeneratorAction::ZIIIPrimaryGeneratorAction(
                                               ZIIIDetectorConstruction* ZIIIDC)
:ZIIIDetector(ZIIIDC)
{
  particleGun = new G4GeneralParticleSource ();
  //  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,20.*cm));
}

ZIIIPrimaryGeneratorAction::~ZIIIPrimaryGeneratorAction()
{
  delete particleGun;
}

void ZIIIPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  // 
  //!  G4double x0 = -0.5*(ZIIIDetector->GetWorldSizeX());
  G4double x0 = 0.*cm, y0 = 0.*cm, z0 = 20.*cm;
  if (rndmFlag == "on")
    {y0 = 10.0*cm*(G4UniformRand()-0.5);
     z0 = 10.0*cm*(G4UniformRand()-0.5);
//!     {y0 = (ZIIIDetector->GetCalorSizeYZ())*(G4UniformRand()-0.5);
//!      z0 = (ZIIIDetector->GetCalorSizeYZ())*(G4UniformRand()-0.5);
     } 
  //  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  particleGun->GeneratePrimaryVertex(anEvent);
}




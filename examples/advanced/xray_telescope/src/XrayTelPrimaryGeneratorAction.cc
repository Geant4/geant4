// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelPrimaryGeneratorAction.cc                *
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

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"

#include "XrayTelPrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelPrimaryGeneratorAction::XrayTelPrimaryGeneratorAction()
{
  particleGun = new G4GeneralParticleSource ();
}

XrayTelPrimaryGeneratorAction::~XrayTelPrimaryGeneratorAction()
{
  delete particleGun;
}

void XrayTelPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->GeneratePrimaryVertex(anEvent);
}




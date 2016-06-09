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
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelPrimaryGeneratorAction.hh                *
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
// - First implementation of X-ray Telescope advanced example.
// - Based on Chandra and XMM models
//
//
// **********************************************************************

#ifndef XrayTelPrimaryGeneratorAction_h
#define XrayTelPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4GeneralParticleSource;
class G4Event;

class XrayTelPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  XrayTelPrimaryGeneratorAction();
  ~XrayTelPrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);

private:
  G4GeneralParticleSource* particleGun;
  
};

#endif




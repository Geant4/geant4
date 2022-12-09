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
/// \file VG01PrimaryGeneratorAction.hh
/// \brief Definition of the VG01PrimaryGeneratorAction class

//  Generator with multiple modes: uniform, free, axis, fixed)
//
//  Primary generator that can be configured using modes to either:
//    - 'uniform' direction secondaries       (default)
//    - 'free'    which is configured by UI - initialised in con/tor
//    - 'axis'    directing secondaries by rotation to x, y & z axis
//    - 'fixed'   along x axis, and near it ( y = 0.1 or z = 0.1)
// 
//  Authors: J. Apostolakis & S. Wenzel (CERN)  2018-2021
//
//  Started from FullCMS code by Mihaly Novak (CERN) 2017  

#ifndef _G01PRIMARYGENERATORACTION_H_
#define _G01PRIMARYGENERATORACTION_H_

#include "G4VUserPrimaryGeneratorAction.hh"

#include "globals.hh"

class G4Event;
class G4ParticleGun;

/// Configurable primary generator action to demonstrate use of 
///   VecGeom Navigation

class VG01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{

  public:
   enum EGeneratorMode { kFixedMode, kUniformMode, kAxisMode,  kFreeMode};
  //                     x Axis       random dir    x, y, z     UI directed

    VG01PrimaryGeneratorAction( EGeneratorMode mode = kUniformMode );
   ~VG01PrimaryGeneratorAction();

   void GeneratePrimaries(G4Event* anEvent) override;

   void SetGeneratorMode(EGeneratorMode val) { fMode = val;} 
   EGeneratorMode GetGeneratorMode() { return fMode;}
  private:

    G4ParticleGun* fParticleGun;
    EGeneratorMode  fMode;
};

#endif

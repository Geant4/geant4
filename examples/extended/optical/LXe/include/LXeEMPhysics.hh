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
// $Id: LXeEMPhysics.hh 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/LXe/include/LXeEMPhysics.hh
/// \brief Definition of the LXeEMPhysics class
//
//
#ifndef LXeEMPhysics_h
#define LXeEMPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

class LXeEMPhysics : public G4VPhysicsConstructor
{
  public:

    LXeEMPhysics(const G4String& name ="EM");
    virtual ~LXeEMPhysics();

  public:

    // This method will be invoked in the Construct() method.
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    virtual void ConstructProcess();

  protected:

   // Gamma physics
    G4PhotoElectricEffect* fPhotoEffect;
    G4ComptonScattering* fComptonEffect;
    G4GammaConversion* fPairProduction;
 
    // Electron physics
    G4eMultipleScattering* fElectronMultipleScattering;
    G4eIonisation* fElectronIonisation;
    G4eBremsstrahlung* fElectronBremsStrahlung;
 
    //Positron physics
    G4eMultipleScattering* fPositronMultipleScattering;
    G4eIonisation* fPositronIonisation;
    G4eBremsstrahlung* fPositronBremsStrahlung;
    G4eplusAnnihilation* fAnnihilation;
};

#endif

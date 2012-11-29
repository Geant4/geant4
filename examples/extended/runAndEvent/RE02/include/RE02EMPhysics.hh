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
/// \file runAndEvent/RE02/include/RE02EMPhysics.hh
/// \brief Definition of the RE02EMPhysics class
//
// $Id$
// --------------------------------------------------------------
//
// 09-Oct-2003 Chhange gamma, electron, positorn process T. Koi

#ifndef RE02EMPhysics_h
#define RE02EMPhysics_h 1

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

//
/// User electromagnetic physics constructor
///
///  applys related processes to gamma and e-/+
///
/// - void ConstructParticle()
///     does nothing
///
/// - void ConstructProcess()
///     adds processes to each particles
///     gamma :
///       G4GammaConversion, G4ComptonScattering and G4PHotoElectricEffect
///     electron : 
///       G4eMultipleScattering, G4eIonisation and G4eBremsstrahlung
///     positron :
///       G4eMultipleScattering, G4eIonisation, G4eBremsstrahlung and
///       G4ePlusAnnihilation
//
class RE02EMPhysics : public G4VPhysicsConstructor
{
  public:
    RE02EMPhysics(const G4String& name ="EM");
    virtual ~RE02EMPhysics();

  public:
    // This method will be invoked in the Construct() method.
    // each particle type will be instantiated
    virtual void ConstructParticle(){;};

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    virtual void ConstructProcess();

  protected:

};

#endif

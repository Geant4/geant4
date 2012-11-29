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
/// \file runAndEvent/RE02/include/RE02IonPhysics.hh
/// \brief Definition of the RE02IonPhysics class
//
// $Id$
// --------------------------------------------------------------
// 05-Jan-2004 Add G4ionIonisation T. Koi
//

#ifndef RE02IonPhysics_h
#define RE02IonPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4HadronElasticProcess.hh"
#include "G4LElastic.hh"

#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"

#include "G4DeuteronInelasticProcess.hh"
#include "G4LEDeuteronInelastic.hh"

#include "G4TritonInelasticProcess.hh"
#include "G4LETritonInelastic.hh"

#include "G4AlphaInelasticProcess.hh"
#include "G4LEAlphaInelastic.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4hMultipleScattering.hh"

//
/// User ion physics constructor
///
///  applys related processes to ions
///
/// - void ConstructParticle()
///     does nothing
///
/// - void ConstructProcess()
///     adds processes to each particle
///     generic ion :
///       G4HadronInelasticProcess with G4TripathiCrossSection,
///       G4IonsShenCrossSection and G4BinaryLightIonReaction,
///       G4hMultipleScattering and G4ionIonisation
///     deuteron :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4DeuteronInelasticProcess with G4LEDeuteronInelastic,
///       G4hMultipleScattering and G4hIonisation
///     triton :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4TritonInelasticProcess with G4LETritonInelastic,
///       G4hMultipleScattering and G4hIonisation
///     alpha :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4AlphaInelasticProcess with G4LEAlphaInelastic,
///       G4hMultipleScattering and G4hIonisation
///     He3 :
///       G4hMultipleScattering and G4hIonisation
//
class RE02IonPhysics : public G4VPhysicsConstructor
{
  public:
    RE02IonPhysics(const G4String& name="ion");
    virtual ~RE02IonPhysics();

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


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
// $Id: G4IonPhysics.hh,v 1.2 2005/12/02 12:40:04 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4IonBuilder
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4IonPhysics_h
#define G4IonPhysics_h 1

#include "globals.hh"

#include "G4VPhysicsConstructor.hh"

class  G4LElastic;
class  G4HadronElasticProcess;
class  G4DeuteronInelasticProcess;
class  G4LEDeuteronInelastic;
class  G4TritonInelasticProcess;
class  G4LETritonInelastic;
class  G4AlphaInelasticProcess;
class  G4LEAlphaInelastic;

class G4IonPhysics : public G4VPhysicsConstructor
{
public:
  G4IonPhysics(const G4String& name="ions");
  virtual ~G4IonPhysics();

  // This method will be invoked in the Construct() method.
  // each particle type will be instantiated
  virtual void ConstructParticle();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();

private:
  // Elastic Process
  G4LElastic*                 theElasticModel;

  // Generic Ion physics
  G4HadronElasticProcess*     theIonElasticProcess;

  // Deuteron physics
  G4HadronElasticProcess*     theDElasticProcess;
  G4DeuteronInelasticProcess* fDeuteronProcess;
  G4LEDeuteronInelastic*      fDeuteronModel;

  // Triton physics
  G4HadronElasticProcess*     theTElasticProcess;
  G4TritonInelasticProcess*   fTritonProcess;
  G4LETritonInelastic*        fTritonModel;

  // Alpha physics
  G4HadronElasticProcess*     theAElasticProcess;
  G4AlphaInelasticProcess*    fAlphaProcess;
  G4LEAlphaInelastic*         fAlphaModel;

  // He3 physics
  G4HadronElasticProcess*     theHe3ElasticProcess;

  G4bool wasActivated;
};

// 2002 by J.P. Wellisch

#endif


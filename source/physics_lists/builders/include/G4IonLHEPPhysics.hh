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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4IonLHEPPhysics
//
// Author:      A.Ribon 16-Oct-2012
//              Copied from the original G4IonPhysics and renamed.
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4IonLHEPPhysics_h
#define G4IonLHEPPhysics_h 1

#include "globals.hh"

#include "G4VPhysicsConstructor.hh"

class  G4DeuteronInelasticProcess;
class  G4LEDeuteronInelastic;
class  G4TritonInelasticProcess;
class  G4LETritonInelastic;
class  G4AlphaInelasticProcess;
class  G4LEAlphaInelastic;

class G4IonLHEPPhysics : public G4VPhysicsConstructor
{
public:
  G4IonLHEPPhysics(G4int verbose =1);
  //obsolete
  G4IonLHEPPhysics(const G4String& name);
  virtual ~G4IonLHEPPhysics();

  // This method will be invoked in the Construct() method.
  // each particle type will be instantiated
  virtual void ConstructParticle();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();

private:

  // Deuteron physics
  G4DeuteronInelasticProcess* fDeuteronProcess;
  G4LEDeuteronInelastic*      fDeuteronModel;

  // Triton physics
  G4TritonInelasticProcess*   fTritonProcess;
  G4LETritonInelastic*        fTritonModel;

  // Alpha physics
  G4AlphaInelasticProcess*    fAlphaProcess;
  G4LEAlphaInelastic*         fAlphaModel;


  G4bool wasActivated;
};

// 2002 by J.P. Wellisch

#endif


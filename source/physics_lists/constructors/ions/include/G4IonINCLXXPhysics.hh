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
//---------------------------------------------------------------------------
//
// ClassName:   G4IonINCLXXBuilder
//
// Author:      D. Mancusi 23.03.2012
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4IonINCLXXPhysics_h
#define G4IonINCLXXPhysics_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"

#include <vector>

class G4INCLXXInterface;
class G4FTFBuilder;
class G4ParticleDefinition;
class G4HadronicInteraction;
class G4VCrossSectionDataSet;

class G4IonINCLXXPhysics : public G4VPhysicsConstructor
{
public:
  G4IonINCLXXPhysics(G4int ver = 0);
  G4IonINCLXXPhysics(const G4String& name, G4int ver = 0);
  virtual ~G4IonINCLXXPhysics();

  // This method will be invoked in the Construct() method.
  // each particle type will be instantiated
  virtual void ConstructParticle();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();

private:

  void AddProcess(const G4String&,
		  G4ParticleDefinition*, 
		  G4HadronicInteraction*,
		  G4HadronicInteraction*,
		  G4VCrossSectionDataSet*);

  static G4ThreadLocal G4INCLXXInterface* theINCLXXDeuteron;
  static G4ThreadLocal G4INCLXXInterface* theINCLXXTriton;
  static G4ThreadLocal G4INCLXXInterface* theINCLXXHe3;
  static G4ThreadLocal G4INCLXXInterface* theINCLXXAlpha;
  static G4ThreadLocal G4INCLXXInterface* theINCLXXIons;
  static G4ThreadLocal G4FTFBuilder* theFTFPBuilder;

  G4double emaxINCLXX;
  G4double deltaE;

  G4int  verbose;
};

#endif


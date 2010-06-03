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
// $Id: G4IonQMDPhysics.hh,v 1.2 2010-06-03 15:03:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4IonBinaryCascadeBuilder
//   Created from G4IonBinaryCascadePhysics
//
// Author:      G.Folger
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4IonQMDPhysics_h
#define G4IonQMDPhysics_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"

#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"

#include <vector>

class G4HadronInelasticProcess;
class G4HadronicInteraction;
class G4BinaryLightIonReaction;
class G4QMDReaction;
class G4TripathiLightCrossSection;
class G4TripathiCrossSection;
class G4IonsShenCrossSection;


class G4IonQMDPhysics : public G4VPhysicsConstructor
{
public:
  G4IonQMDPhysics(G4int verb = 0);
  G4IonQMDPhysics(const G4String& name, G4int ver = 0);
  virtual ~G4IonQMDPhysics();

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
		  G4BinaryLightIonReaction*,
		  G4QMDReaction*,
		  G4HadronicInteraction*);

  std::vector<G4HadronInelasticProcess*> p_list;
  std::vector<G4HadronicInteraction*> model_list;

  G4TripathiCrossSection* fTripathi;
  G4TripathiLightCrossSection* fTripathiLight;
  G4IonsShenCrossSection* fShen;

  G4LEDeuteronInelastic*  fLEDModel;
  G4LETritonInelastic*    fLETModel;
  G4LEAlphaInelastic*     fLEAModel;

  G4double eminBIC;
 
  G4double eminQMD;
  G4double emaxQMD;

  G4double emaxLHEP;

  G4double overlap;
   
  G4int  verbose;
  G4bool wasActivated;
};


#endif


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
// $Id: G4IonBinaryCascadePhysics.hh,v 1.2 2006-06-06 16:47:44 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4IonBinaryCascadeBuilder
//
// Author:      V.Ivanchenko 09.04.2006
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4IonBinaryCascadePhysics_h
#define G4IonBinaryCascadePhysics_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"
#include <vector>

class G4HadronInelasticProcess;
class G4HadronicInteraction;
class G4TripathiLightCrossSection;
class G4TripathiCrossSection;
class G4IonsShenCrossSection;

class G4IonBinaryCascadePhysics : public G4VPhysicsConstructor
{
public:
  G4IonBinaryCascadePhysics(const G4String& name="ions", G4int verb = 0);
  virtual ~G4IonBinaryCascadePhysics();

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
		  G4HadronicInteraction*);

  std::vector<G4HadronInelasticProcess*> p_list;
  std::vector<G4HadronicInteraction*> model_list;

  G4TripathiCrossSection* fTripathi;
  G4TripathiLightCrossSection* fTripathiLight;
  G4IonsShenCrossSection* fShen;

  G4double emax;
  G4double emaxLHEP;
  G4double eminBIC;

  G4int  verbose;
  G4bool wasActivated;
};


#endif


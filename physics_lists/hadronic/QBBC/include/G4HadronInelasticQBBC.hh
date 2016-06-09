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
// $Id: G4HadronInelasticQBBC.hh,v 1.5 2006/07/05 16:12:43 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-01-patch-01 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronInelasticQBBC
//
// Author: 11 April 2006 V. Ivanchenko
//
// Modified:
// 05.07.2006 V.Ivanchenko fix problem of initialisation of HP
//
//----------------------------------------------------------------------------
//

#ifndef G4HadronInelasticQBBC_h
#define G4HadronInelasticQBBC_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4FTFModel.hh"
#include "G4CascadeInterface.hh"

#include "G4PiNuclearCrossSection.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPFissionData.hh"

class G4HadronicProcess;
class G4TheoFSGenerator;
class G4StringChipsParticleLevelInterface;
class G4ExcitedStringDecay;
class G4HadronProcessStore;

class G4HadronInelasticQBBC : public G4VPhysicsConstructor
{
public: 

  G4HadronInelasticQBBC(const G4String& name = "inelastic",
			G4int ver = 1, G4bool ftf = false, G4bool bert = false,
			G4bool chips = false, G4bool hp = false);

  virtual ~G4HadronInelasticQBBC();

public: 

  // This method will be invoked in the Construct() method. 
  // each particle type will be instantiated
  virtual void ConstructParticle();
 
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  virtual void ConstructProcess();

private:

  void Register(G4ParticleDefinition*, G4HadronicProcess*, 
		G4HadronicInteraction*, const G4String&);

  G4PiNuclearCrossSection thePiCross;
  G4ProtonInelasticCrossSection  theXSecP;
  G4NeutronInelasticCrossSection theXSecN;

  G4NeutronHPInelasticData*  theHPXSecI;
  G4NeutronHPCaptureData*    theHPXSecC;
  G4NeutronHPFissionData*    theHPXSecF;

  G4HadronProcessStore* store;

  G4StringChipsParticleLevelInterface * theCHIPSCascade;
  G4QGSModel< G4QGSParticipants > * theQGStringModel;
  G4ExcitedStringDecay* theQGStringDecay;
  G4ExcitedStringDecay* theFTFStringDecay;
  G4FTFModel*           theFTFStringModel;

  G4int    verbose;
  G4bool   ftfFlag;
  G4bool   bertFlag;
  G4bool   chipsFlag;
  G4bool   hpFlag;
  G4bool   wasActivated;
};

#endif

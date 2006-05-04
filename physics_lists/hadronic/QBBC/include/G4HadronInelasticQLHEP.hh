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
// $Id: G4HadronInelasticQLHEP.hh,v 1.1 2006-05-04 16:48:39 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronInelasticQLHEP
//
// Author: 11 April 2006 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4HadronInelasticQLHEP_h
#define G4HadronInelasticQLHEP_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4FTFModel.hh"

#include "G4PiNuclearCrossSection.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPFissionData.hh"

#include <vector>

class G4HadronicProcess;
class G4TheoFSGenerator;
class G4GeneratorPrecompoundInterface;
class G4PreCompoundModel;
class G4ExcitedStringDecay;
class G4HadronProcessStore;

class G4HadronInelasticQLHEP : public G4VPhysicsConstructor
{
public: 

  G4HadronInelasticQLHEP(const G4String& name = "inelastic",
			G4int ver = 1, G4bool qgs = false, G4bool bert = false,
			G4bool bic = false, G4bool hp = false);

  virtual ~G4HadronInelasticQLHEP();

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

  void AddLHEP(G4ParticleDefinition*,
	       G4HadronicProcess*,
	       G4double emin, G4double emax);

  G4PiNuclearCrossSection thePiCross;
  G4ProtonInelasticCrossSection  theXSecP;
  G4NeutronInelasticCrossSection theXSecN;

  G4NeutronHPInelasticData  theHPXSecI;
  G4NeutronHPCaptureData    theHPXSecC;
  G4NeutronHPFissionData    theHPXSecF;

  G4HadronProcessStore* store;

  G4GeneratorPrecompoundInterface * theCascade;
  G4PreCompoundModel * thePreEquilib;
  G4QGSModel< G4QGSParticipants > * theQGStringModel;
  G4ExcitedStringDecay* theQGStringDecay;

  G4int    verbose;
  G4bool   qgsFlag;
  G4bool   bertFlag;
  G4bool   bicFlag;
  G4bool   hpFlag;
  G4bool   wasActivated;
};

#endif

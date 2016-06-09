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
//    **************************************
//    *                                    *
//    *    RemSimHadronicPhysics.hh        *
//    *                                    *
//    **************************************
//
// $Id: RemSimHadronicPhysics.hh,v 1.2 2004/05/22 12:57:04 guatelli Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// Author : Susanna Guatelli, guatelli@ge.infn.it
// 

#ifndef RemSimHadronicPhysics_h
#define RemSimHadronicPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include "G4QGSParticipants.hh"
#include "G4QGSModel.hh"

class G4ProtonInelasticCrossSection;
class G4NeutronInelasticCrossSection;
class G4PiNuclearCrossSection;
class G4TripathiCrossSection;
class G4IonsShenCrossSection;

class G4LElastic;
class G4CascadeInterface;
class G4BinaryLightIonReaction;
class G4LEProtonInelastic;
class G4LENeutronInelastic;
class G4LEPionPlusInelastic;
class G4LEPionMinusInelastic;
class G4LEAlphaInelastic;
class G4LFission;
class G4LCapture;

class G4TheoFSGenerator;
class G4GeneratorPrecompoundInterface;
class G4ExcitationHandler;
class G4PreCompoundModel;
class G4QGSMFragmentation;
class G4ExcitedStringDecay;

class RemSimHadronicPhysics: public G4VPhysicsConstructor 
{
  public:
    RemSimHadronicPhysics(const G4String& name = "hadronic");
    virtual ~RemSimHadronicPhysics();

  protected:
    // Construct particle and physics
    void ConstructParticle(){};
    void ConstructProcess();

  private:

    G4ProtonInelasticCrossSection* proton_XC;
    G4NeutronInelasticCrossSection* neutron_XC;
    G4PiNuclearCrossSection* pion_XC;
 
    G4TripathiCrossSection* tripathi;
    G4IonsShenCrossSection* shen;

    G4LElastic* elastic_model;
    G4CascadeInterface* bertini_model;
    G4BinaryLightIonReaction* binary_ion_model;
    
    G4LEProtonInelastic* LEP_proton_model;
    G4LENeutronInelastic* LEP_neutron_model;
    G4LEPionPlusInelastic* LEP_pip_model;
    G4LEPionMinusInelastic* LEP_pim_model;
    G4LEAlphaInelastic* LEP_alpha_model;
    G4LFission* nfission_model;
    G4LCapture* ncapture_model;

    G4TheoFSGenerator* QGSP_model;
    G4GeneratorPrecompoundInterface* theCascade;
    G4ExcitationHandler* theHandler;
    G4PreCompoundModel* thePreEquilib;
    G4QGSMFragmentation* theFragmentation;
    G4ExcitedStringDecay* theStringDecay;
    G4QGSModel<G4QGSParticipants> theStringModel;
};
#endif









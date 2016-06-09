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
//    **************************************
//    *                                    *
//    *    RemSimHadronicBinary.hh        *
//    *                                    *
//    **************************************
//
// $Id: RemSimHadronicBinary.hh,v 1.4 2006/06/29 16:22:47 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// Author : Susanna Guatelli, guatelli@ge.infn.it
// 

#ifndef RemSimHadronicBinary_h
#define RemSimHadronicBinary_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"
#include "G4ExcitationHandler.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSModel.hh"

class G4ProtonInelasticCrossSection;
class G4NeutronInelasticCrossSection;
class G4PiNuclearCrossSection;
class G4TripathiCrossSection;
class G4IonsShenCrossSection;
class G4BinaryCascade;
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

class RemSimHadronicBinary: public G4VPhysicsConstructor 
{
  public:
    RemSimHadronicBinary(const G4String& name = "hadronic-binary");
    virtual ~RemSimHadronicBinary();

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
    G4BinaryLightIonReaction* binary_ion_model;
    G4BinaryCascade* binarycascade_model;
    G4LEProtonInelastic* LEP_proton_model;
    G4LENeutronInelastic* LEP_neutron_model;
    G4LEPionPlusInelastic* LEP_pip_model;
    G4LEPionMinusInelastic* LEP_pim_model;
    G4LEAlphaInelastic* LEP_alpha_model;
    G4LFission* nfission_model;
    G4LCapture* ncapture_model;

    G4TheoFSGenerator* QGSP_model;
    G4GeneratorPrecompoundInterface* theCascade;
    G4ExcitationHandler theHandler;
    G4PreCompoundModel* thePreEquilib;
    G4QGSMFragmentation* theFragmentation;
    G4ExcitedStringDecay* theStringDecay;
    G4QGSModel<G4QGSParticipants> theStringModel;
};
#endif









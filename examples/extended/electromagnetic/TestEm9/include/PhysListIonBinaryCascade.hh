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
// $Id: PhysListIonBinaryCascade.hh,v 1.2 2003/12/05 11:22:45 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// ------------------------------------------------------------
//	GEANT 4 class header file
// Class Description:
//      This class is an derived class of G4VPhysicsConstructor
//      It is provide PhysicsList for Binary Cascade for
//      protons and neutrons with the energy E<3 GeV
//
// ------------------------------------------------------------
//	History
//        Created:       18.11.03  V.Ivanchenko
//
//        Modified:
//
// ------------------------------------------------------------
//
#ifndef PhysListIonBinaryCascade_h
#define PhysListIonBinaryCascade_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

class G4HadronInelasticProcess;

class PhysListIonBinaryCascade : public G4VPhysicsConstructor
{
  public:
    PhysListIonBinaryCascade(const G4String& name = "binary_ion");
    virtual ~PhysListIonBinaryCascade();

  public:
    // This method will be invoked in the Construct() method.
    // each particle type will be instantiated
    void ConstructParticle() {};

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    void ConstructProcess();

  private:

    G4DeuteronInelasticProcess      theIPdeuteron;
    G4TritonInelasticProcess        theIPtriton;
    G4AlphaInelasticProcess         theIPalpha;
    G4HadronInelasticProcess*       theIPHe3;
    G4HadronInelasticProcess*       theIPIonC12;
    G4HadronInelasticProcess*       theIPGenericIon;

};

#endif
